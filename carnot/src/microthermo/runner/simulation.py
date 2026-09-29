from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any
import copy
import json
import math
from pathlib import Path
import numpy as np

from ..core.boundaries import BoundaryKind, SegmentWall, Wall, resolve_wall
from ..core.apparatus import ApparatusComponent
from ..core.ccd import body_pair_toi
from ..core.broadphase import all_pairs, swept_pairs
from ..core.numeric import (fixed_wall_toi, numba_available, polygon_wall_gap,
                            reach_mask, wall_lower_bounds)
from ..core.geometry import (closest_point_segment, convex_separation,
                             disc_polygon_separation, polygon_segment_witnesses,
                             world_polygon)
from ..core.contacts import (apply_impulse, contact_velocity, inverse_effective_mass,
                             elastic_cluster_impulses, resolve_elastic, resolve_energy_step)
from ..core.events import Contact, InteractionRecord, TOIResult, TOIStatus
from ..core.state import BodyArrays, Shape, Tolerances
from ..measurements.ledger import EnergyLedger
from ..measurements.cycles import CycleMarker
from ..measurements.observables import rotational_temperature, translational_temperature


@dataclass
class Portal:
    point: np.ndarray
    normal: np.ndarray
    half_length: float
    host_id: int = -1
    inside_label: int = 1
    outside_label: int = 0
    delta_u: float = 0.0
    inside_energy: float = 0.0


@dataclass
class World:
    bodies: BodyArrays
    walls: list[Wall]
    portals: list[Portal] = field(default_factory=list)
    mechanism: Any = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def energy(self) -> float:
        value = self.bodies.kinetic_energy()
        if self.mechanism is not None:
            value += self.mechanism.energy()
        return value


@dataclass(frozen=True)
class Snapshot:
    time: float
    position: np.ndarray
    velocity: np.ndarray
    angle: np.ndarray
    omega: np.ndarray
    energy: float
    translational_temperature: float
    rotational_temperature: float
    branch: str | None
    occupancy: tuple[int, ...]
    event_count: int
    energy_residual: float
    apparatus: tuple[ApparatusComponent, ...] = ()
    shaft_phase: float | None = None


@dataclass
class Checkpoint:
    time: float
    bodies: BodyArrays
    mechanism: Any
    rng_state: dict
    ledger: EnergyLedger
    memberships: np.ndarray
    event_count: int
    pending: tuple[Checkpoint, TOIResult, float] | None = None
    cluster_count: int = 0
    max_cluster_residual: float = 0.0
    history: tuple[InteractionRecord, ...] | None = None
    cycle_markers: tuple[CycleMarker, ...] = ()


class NumericalFailure(RuntimeError):
    pass


class Simulation:
    """Deterministic bounded-horizon event-driven reference scheduler."""

    def __init__(self, world: World, seed: int = 123, max_horizon: float = 0.05,
                 tolerances: Tolerances = Tolerances(), pair_search: str = "grid",
                 numeric_backend: str = "auto", wall_search: str = "bounded",
                 wall_kernel: str = "auto", penetration_kernel: str = "auto"):
        if pair_search not in ("grid", "all"):
            raise ValueError("pair_search must be grid or all")
        if numeric_backend not in ("auto", "python", "numba"):
            raise ValueError("numeric_backend must be auto, python, or numba")
        if numeric_backend == "numba" and not numba_available():
            raise RuntimeError("Numba backend requires installation of microthermo[accel]")
        if wall_search not in ("bounded", "all"):
            raise ValueError("wall_search must be bounded or all")
        if wall_kernel not in ("auto", "python"):
            raise ValueError("wall_kernel must be auto or python")
        if penetration_kernel not in ("auto", "python"):
            raise ValueError("penetration_kernel must be auto or python")
        self.pair_search = pair_search
        self.numeric_backend = numeric_backend
        self.wall_search = wall_search
        self.wall_kernel = wall_kernel
        self.penetration_kernel = penetration_kernel
        self.world = world
        self.time = 0.0
        self.rng = np.random.default_rng(seed)
        self.seed = seed
        self.max_horizon = max_horizon
        self.tol = tolerances
        self.events: list[InteractionRecord] = []
        self.samples: list[Snapshot] = []
        self.memberships = np.zeros(world.bodies.n, dtype=np.int16)
        self.ledger = EnergyLedger(self.energy())
        self.event_count = 0
        self.max_penetration = 0.0
        self.ccd_refinements = 0
        self.ccd_failures = 0
        self.cluster_count = 0
        self.max_cluster_residual = 0.0
        self.failure_checkpoint: Checkpoint | None = None
        self.failure_diagnostic: dict[str, Any] | None = None
        self.cycle_markers: list[CycleMarker] = []
        shaft=getattr(world.mechanism,"shaft",None)
        self._controlled_phase_origin=(shaft.phi if shaft is not None and
                                       shaft.prescribed_omega is not None else None)
        # An interval is searched from a fixed physical origin. Intermediate
        # snapshots may inspect it without changing the next collision query.
        self._pending: tuple[Checkpoint, TOIResult, float] | None = None
        gap, participants = self._penetration()
        if gap < -self.tol.geometry:
            raise ValueError(f"invalid overlapping initial state: {participants}, gap={gap}")
        mechanism=self.world.mechanism
        if hasattr(mechanism,"cam") and hasattr(mechanism,"shaft"):
            phase=mechanism.shaft.phi
            if abs((phase+math.pi)%(2*math.pi)-math.pi) <= 1e-10:
                self.cycle_markers.append(self._cycle_marker())

    def energy(self) -> float:
        value=self.world.energy()
        levels={p.inside_label:p.inside_energy for p in self.world.portals}
        value += sum(levels.get(int(region),0.) for region in self.memberships)
        return value

    def _cycle_marker(self) -> CycleMarker:
        mechanism=self.world.mechanism
        ledger=self.ledger
        return CycleMarker(
            time=self.time,phase=mechanism.shaft.phi,
            total_energy=self.energy(),
            gas_energy=self.world.bodies.kinetic_energy(),
            apparatus_energy=mechanism.energy(),
            heat_hot=ledger.heat_hot.value,heat_cold=ledger.heat_cold.value,
            heat_other=ledger.heat_other.value,work_on=ledger.work_on.value,
            piston_work_on_gas=ledger.piston_work_on_gas.value,
            motor_work=mechanism.motor_work,load_output=ledger.load_output.value,
            event_count=self.event_count)

    def checkpoint(self, include_history: bool = True) -> Checkpoint:
        return Checkpoint(self.time, self.world.bodies.copy(), copy.deepcopy(self.world.mechanism),
                          copy.deepcopy(self.rng.bit_generator.state), copy.deepcopy(self.ledger),
                          self.memberships.copy(), self.event_count,
                          copy.deepcopy(self._pending),self.cluster_count,
                          self.max_cluster_residual,
                          tuple(copy.deepcopy(self.events)) if include_history else None,
                          tuple(self.cycle_markers))

    def restore(self, cp: Checkpoint) -> None:
        self.time, self.world.bodies, self.world.mechanism = cp.time, cp.bodies.copy(), copy.deepcopy(cp.mechanism)
        self.rng.bit_generator.state = copy.deepcopy(cp.rng_state)
        self.ledger, self.memberships, self.event_count = copy.deepcopy(cp.ledger), cp.memberships.copy(), cp.event_count
        self.events = (list(copy.deepcopy(cp.history)) if cp.history is not None
                       else self.events[:cp.event_count])
        self._pending = copy.deepcopy(cp.pending)
        self.cluster_count = cp.cluster_count
        self.max_cluster_residual = cp.max_cluster_residual
        self.cycle_markers = list(cp.cycle_markers)
        if self.world.mechanism is not None:
            for wall in self.world.walls:
                if hasattr(wall, "mechanism"):
                    wall.mechanism = self.world.mechanism
            clock = self.world.metadata.get("clock")
            if clock is not None:
                clock[0] = self.time

    def _restore_interval_origin(self, cp: Checkpoint) -> None:
        self.restore(cp)
        self._pending = None

    def _wall_toi(self, body: int, wi: int, horizon: float,
                  wall_speed_bound: float | None = None,
                  wall_start: np.ndarray | None = None) -> TOIResult:
        wall, s = self.world.walls[wi], self.world.bodies
        if (wall_start is not None and not isinstance(wall, SegmentWall)):
            inward = wall.inward_normal
            gap_bound = float((s.pos[body] - wall_start) @ inward - s.radius[body])
            body_speed = abs(float(s.vel[body] @ inward))
            if gap_bound > horizon * (body_speed + wall_speed_bound) + self.tol.geometry:
                return TOIResult(TOIStatus.NO_COLLISION)
        if isinstance(wall, SegmentWall):
            if s.shape[body] != Shape.DISC:
                return self._polygon_segment_wall_toi(body, wi, horizon)
            return self._segment_wall_toi(body, wi, horizon)
        if s.shape[body] != Shape.DISC:
            return self._polygon_wall_toi(body, wi, horizon, wall_speed_bound)
        inward = wall.inward_normal
        t = 0.0
        bound = (np.linalg.norm(s.vel[body]) +
                 (wall_speed_bound if wall_speed_bound is not None else
                  wall.speed_bound(self.time, self.time+horizon)))
        if bound <= self.tol.velocity:
            return TOIResult(TOIStatus.NO_COLLISION)
        local_time_tol = min(self.tol.time, self.tol.geometry/bound)
        previous = 0.0
        for _ in range(1024):
            center = s.pos[body] + s.vel[body]*t
            wp = wall.point_at(self.time+t)
            gap = float((center-wp) @ inward - s.radius[body])
            rel = float((s.vel[body]-wall.velocity_at(self.time+t)) @ inward)
            if gap <= self.tol.geometry:
                if t == 0.0 and rel >= -self.tol.velocity:
                    return TOIResult(TOIStatus.NO_COLLISION, reason="touching and separating")
                lo, hi = previous, t
                for _ in range(60):
                    if hi-lo <= local_time_tol:
                        break
                    mid = 0.5*(lo+hi)
                    c = s.pos[body]+s.vel[body]*mid
                    g = float((c-wall.point_at(self.time+mid)) @ inward-s.radius[body])
                    if g <= self.tol.geometry: hi = mid
                    else: lo = mid
                center = s.pos[body]+s.vel[body]*hi
                point = center-s.radius[body]*inward
                return TOIResult(TOIStatus.COLLISION, hi, hi-lo,
                                 Contact(body, None, tuple(point), tuple(-inward), boundary=wi))
            previous = t
            step = 0.8*gap/bound
            if step <= local_time_tol:
                return TOIResult(TOIStatus.INDETERMINATE, reason="wall advancement stalled")
            t += step
            if t > horizon:
                return TOIResult(TOIStatus.NO_COLLISION)
        return TOIResult(TOIStatus.INDETERMINATE, reason="wall iteration limit")

    def _polygon_wall_toi(self, body: int, wi: int, horizon: float,
                          wall_speed_bound: float | None = None) -> TOIResult:
        wall, s = self.world.walls[wi], self.world.bodies
        inward = wall.inward_normal
        speed = (np.linalg.norm(s.vel[body]) + abs(s.omega[body])*s.radius[body]
                 + (wall_speed_bound if wall_speed_bound is not None else
                    wall.speed_bound(self.time, self.time+horizon)))
        if (self.wall_kernel == "auto" and
            type(wall).point_at is Wall.point_at and
            type(wall).velocity_at is Wall.velocity_at and
            (self.numeric_backend == "numba" or
             (self.numeric_backend == "auto" and numba_available()))):
            status, time, error, px, py, feature, reason = fixed_wall_toi(
                s, body, wall, wall.point_at(self.time), speed, horizon, self.tol)
            if status == 1:
                return TOIResult(TOIStatus.COLLISION, time, error,
                                 Contact(body, None, (px, py), tuple(-inward),
                                         feature_a=2*feature, boundary=wi))
            if status == 0:
                return TOIResult(TOIStatus.NO_COLLISION)
            reasons = ("", "initial polygon-wall overlap",
                       "touching polygon-wall contact cannot be isolated",
                       "new polygon-wall contact during separation",
                       "polygon-wall contact did not separate",
                       "polygon-wall advancement stalled",
                       "polygon-wall iteration limit")
            return TOIResult(TOIStatus.INDETERMINATE, reason=reasons[reason])
        if speed <= self.tol.velocity:
            return TOIResult(TOIStatus.NO_COLLISION)
        local_time_tol=min(self.tol.time,self.tol.geometry/speed)

        def gap_at(t):
            vertices = world_polygon(s.polygons[body], s.pos[body]+s.vel[body]*t,
                                     s.angle[body]+s.omega[body]*t)
            projection = vertices @ inward
            vertex = int(np.argmin(projection))
            gap = float(projection[vertex]-wall.point_at(self.time+t) @ inward)
            return gap, vertices[vertex], vertex

        t, previous = 0.0, 0.0
        for _ in range(1024):
            gap, point, feature = gap_at(t)
            if gap <= self.tol.geometry:
                if gap < -self.tol.geometry and t == 0:
                    return TOIResult(TOIStatus.INDETERMINATE, reason="initial polygon-wall overlap")
                relative = (s.vel[body]-wall.velocity_at(self.time+t)
                            + s.omega[body]*np.array([-(point-s.pos[body]-s.vel[body]*t)[1],
                                                       (point-s.pos[body]-s.vel[body]*t)[0]]))
                if t <= local_time_tol and relative @ inward >= -self.tol.velocity:
                    t=min(horizon,max(16*local_time_tol,16*self.tol.geometry/speed))
                    if t <= local_time_tol:
                        return TOIResult(TOIStatus.INDETERMINATE, reason="touching polygon-wall contact cannot be isolated")
                    next_gap=gap_at(t)[0]
                    if next_gap < -self.tol.geometry:
                        return TOIResult(TOIStatus.INDETERMINATE, reason="new polygon-wall contact during separation")
                    if next_gap <= self.tol.geometry:
                        return TOIResult(TOIStatus.INDETERMINATE, reason="polygon-wall contact did not separate")
                    previous=t
                    continue
                lo, hi = previous, t
                for _ in range(60):
                    if hi-lo <= local_time_tol: break
                    mid = .5*(lo+hi)
                    if gap_at(mid)[0] <= self.tol.geometry: hi=mid
                    else: lo=mid
                _, point, feature = gap_at(hi)
                return TOIResult(TOIStatus.COLLISION, hi, hi-lo,
                                 Contact(body, None, tuple(point), tuple(-inward),
                                         feature_a=2*feature, boundary=wi))
            previous=t
            step=.8*gap/speed
            if step <= local_time_tol:
                return TOIResult(TOIStatus.INDETERMINATE, reason="polygon-wall advancement stalled")
            t+=step
            if t>horizon: return TOIResult(TOIStatus.NO_COLLISION)
        return TOIResult(TOIStatus.INDETERMINATE, reason="polygon-wall iteration limit")

    def _segment_wall_toi(self, body: int, wi: int, horizon: float) -> TOIResult:
        wall = self.world.walls[wi]
        s = self.world.bodies
        p, v = s.pos[body], s.vel[body]
        a, b = wall.start, wall.end
        edge = b-a; length = np.linalg.norm(edge); tangent = edge/length
        normal = np.array([-tangent[1], tangent[0]])
        effective_r = s.radius[body]+wall.thickness
        candidates: list[tuple[float, np.ndarray, np.ndarray]] = []
        # Hits on either face of the finite segment.
        signed = float((p-a)@normal); vn = float(v@normal)
        if abs(vn) > self.tol.velocity:
            for target in (-effective_r, effective_r):
                t = (target-signed)/vn
                if -self.tol.time <= t <= horizon+self.tol.time:
                    center = p+max(0.,t)*v
                    along = float((center-a)@tangent)
                    if -self.tol.geometry <= along <= length+self.tol.geometry:
                        closest = a+min(length,max(0.,along))*tangent
                        d = closest-center; norm=np.linalg.norm(d)
                        if norm > 0:
                            contact_normal=d/norm
                            if t > self.tol.time or float(v@contact_normal) > self.tol.velocity:
                                candidates.append((max(0.,t), closest, contact_normal))
        # Rounded material endpoints.
        for endpoint in (a,b):
            dp=p-endpoint; qa=float(v@v); qb=float(dp@v); qc=float(dp@dp-effective_r**2)
            if qc <= self.tol.geometry:
                if qb < -self.tol.velocity:
                    d=endpoint-p; candidates.append((0.,endpoint,d/np.linalg.norm(d)))
            elif qa > self.tol.velocity**2 and qb < 0:
                disc=qb*qb-qa*qc
                if disc >= 0:
                    t=qc/(-qb+math.sqrt(max(0.,disc)))
                    if t <= horizon+self.tol.time:
                        center=p+t*v; d=endpoint-center
                        candidates.append((t,endpoint,d/np.linalg.norm(d)))
        if not candidates:
            return TOIResult(TOIStatus.NO_COLLISION)
        t, surface, n = min(candidates,key=lambda x:x[0])
        contact_point = p+v*t+s.radius[body]*n
        return TOIResult(TOIStatus.COLLISION,t,self.tol.time,
                         Contact(body,None,tuple(contact_point),tuple(n),boundary=wi))

    def _polygon_segment_wall_toi(self, body: int, wi: int, horizon: float) -> TOIResult:
        wall, s = self.world.walls[wi], self.world.bodies
        speed=float(np.linalg.norm(s.vel[body])+abs(s.omega[body])*s.radius[body])
        if speed <= self.tol.velocity:
            return TOIResult(TOIStatus.NO_COLLISION)
        local_time_tol=min(self.tol.time,self.tol.geometry/speed)

        def query(t):
            center=s.pos[body]+s.vel[body]*t
            vertices=world_polygon(s.polygons[body],center,s.angle[body]+s.omega[body]*t)
            distance,normal,wp,ws,fp,fs=polygon_segment_witnesses(
                vertices,wall.start,wall.end)
            return distance-wall.thickness,normal,wp,ws,fp,fs,center

        t=previous=0.0
        for _ in range(1024):
            gap,normal,wp,ws,fp,fs,center=query(t)
            if gap <= self.tol.geometry:
                if t == 0.0 and gap < -self.tol.geometry:
                    return TOIResult(TOIStatus.INDETERMINATE,
                                     reason="initial polygon-segment overlap")
                arm=wp-center
                surface_velocity=s.vel[body]+s.omega[body]*np.array([-arm[1],arm[0]])
                closing=float(surface_velocity@normal)
                if closing <= self.tol.velocity:
                    next_t=min(horizon,t+max(16*local_time_tol,
                                             16*self.tol.geometry/speed))
                    if next_t <= t+local_time_tol:
                        return TOIResult(TOIStatus.NO_COLLISION)
                    next_gap=query(next_t)[0]
                    if next_gap < -self.tol.geometry:
                        return TOIResult(TOIStatus.INDETERMINATE,
                                         reason="polygon-segment contact during separation")
                    if next_gap <= self.tol.geometry:
                        return TOIResult(TOIStatus.INDETERMINATE,
                                         reason="polygon-segment grazing contact unresolved")
                    previous=t=next_t
                    continue
                lo,hi=previous,t
                for _ in range(60):
                    if hi-lo <= local_time_tol:
                        break
                    mid=.5*(lo+hi)
                    if query(mid)[0] <= self.tol.geometry:
                        hi=mid
                    else:
                        lo=mid
                _,normal,wp,ws,fp,fs,_=query(hi)
                contact_point=.5*(wp+ws-wall.thickness*normal)
                return TOIResult(TOIStatus.COLLISION,hi,hi-lo,
                                 Contact(body,None,tuple(contact_point),tuple(normal),
                                         feature_a=fp,feature_b=fs,boundary=wi))
            previous=t
            step=.8*gap/speed
            if step <= local_time_tol:
                return TOIResult(TOIStatus.INDETERMINATE,
                                 reason="polygon-segment advancement stalled")
            t+=step
            if t>horizon:
                return TOIResult(TOIStatus.NO_COLLISION)
        return TOIResult(TOIStatus.INDETERMINATE,
                         reason="polygon-segment iteration limit")

    def _portal_toi(self, body: int, pi: int, horizon: float) -> TOIResult:
        portal=self.world.portals[pi]; s=self.world.bodies
        n=portal.normal/np.linalg.norm(portal.normal); tangent=np.array([-n[1],n[0]])
        signed=float((s.pos[body]-portal.point)@n); vn=float(s.vel[body]@n)
        if abs(signed) <= self.tol.geometry or abs(vn) <= self.tol.velocity:
            return TOIResult(TOIStatus.NO_COLLISION)
        t=-signed/vn
        if t < -self.tol.time or t > horizon+self.tol.time:
            return TOIResult(TOIStatus.NO_COLLISION)
        point=s.pos[body]+max(0.,t)*s.vel[body]
        if abs(float((point-portal.point)@tangent)) > portal.half_length:
            return TOIResult(TOIStatus.NO_COLLISION)
        # Orient from the particle toward its direction of crossing.
        cn=-n if vn < 0 else n
        return TOIResult(TOIStatus.COLLISION,max(0.,t),self.tol.time,
                         Contact(body,None,tuple(point),tuple(cn),boundary=-(pi+1)))

    def _earliest(self, horizon: float) -> TOIResult:
        candidates: list[TOIResult] = []
        s = self.world.bodies
        pairs = (swept_pairs(s.pos, s.vel, s.radius, horizon)
                 if self.pair_search == "grid" and s.n >= 32 else all_pairs(s.n))
        if pairs:
            indexed = np.asarray(pairs, dtype=np.int64)
            mask = reach_mask(s.pos, s.vel, s.radius, s.omega, indexed, horizon,
                              self.numeric_backend)
            for (a, b), possible in zip(pairs, mask):
                if not possible:
                    continue
                q = body_pair_toi(s, a, b, horizon, self.tol, self.numeric_backend)
                if q.status == TOIStatus.INDETERMINATE:
                    return TOIResult(q.status, q.time, q.error, q.contact,
                                     f"pair ({a}, {b}): {q.reason}")
                if q.status == TOIStatus.COLLISION:
                    candidates.append(q)
        wall_bounds = ([wall.speed_bound(self.time, self.time+horizon)
                        for wall in self.world.walls]
                       if self.wall_search == "bounded" else None)
        wall_starts = ([wall.point_at(self.time) for wall in self.world.walls]
                       if self.wall_search == "bounded" else None)
        for body in range(s.n):
            for wi in range(len(self.world.walls)):
                q = self._wall_toi(body, wi, horizon,
                                   None if wall_bounds is None else wall_bounds[wi],
                                   None if wall_starts is None else wall_starts[wi])
                if q.status == TOIStatus.INDETERMINATE:
                    return TOIResult(q.status, q.time, q.error, q.contact,
                                     f"wall ({body}, {wi}): {q.reason}")
                if q.status == TOIStatus.COLLISION:
                    candidates.append(q)
            for pi in range(len(self.world.portals)):
                q=self._portal_toi(body,pi,horizon)
                if q.status == TOIStatus.COLLISION:
                    candidates.append(q)
        if not candidates:
            return TOIResult(TOIStatus.NO_COLLISION)
        candidates.sort(key=lambda q: q.time)
        earliest = candidates[0].time
        simultaneous = tuple(q.contact for q in candidates
                             if q.time-earliest <= self.tol.time)
        return TOIResult(TOIStatus.COLLISION, earliest,
                         max(q.error for q in candidates if q.time-earliest <= self.tol.time),
                         simultaneous[0], contacts=simultaneous)

    def _penetration(self) -> tuple[float, tuple[Any, ...] | None]:
        """Return deepest overlap of dynamic bodies and infinite-line walls."""
        s = self.world.bodies
        deepest, participants = 0.0, None
        pairs = (swept_pairs(s.pos, s.vel, s.radius, 0.0)
                 if self.pair_search == "grid" and s.n >= 32 else all_pairs(s.n))
        by_a: list[list[int]] = [[] for _ in range(s.n)]
        for a, b in pairs:
            by_a[a].append(b)
        wall_points = ([wall.point_at(self.time) for wall in self.world.walls]
                       if self.wall_search == "bounded" else None)
        compiled_walls = (self.penetration_kernel == "auto" and wall_points is not None and
                          (self.numeric_backend == "numba" or
                           (self.numeric_backend == "auto" and numba_available())))
        lower_bounds = (wall_lower_bounds(
            s.pos, s.radius, np.asarray(wall_points),
            np.asarray([wall.inward_normal for wall in self.world.walls]))
            if compiled_walls and self.world.walls else None)
        for a in range(s.n):
            for b in by_a[a]:
                if np.linalg.norm(s.pos[b]-s.pos[a]) > s.radius[a]+s.radius[b]:
                    continue
                if s.shape[a] == Shape.DISC and s.shape[b] == Shape.DISC:
                    gap = float(np.linalg.norm(s.pos[b]-s.pos[a])-s.radius[a]-s.radius[b])
                elif s.shape[a] != Shape.DISC and s.shape[b] != Shape.DISC:
                    pa=world_polygon(s.polygons[a],s.pos[a],s.angle[a])
                    pb=world_polygon(s.polygons[b],s.pos[b],s.angle[b])
                    gap=convex_separation(pa,pb)[0]
                else:
                    disc, poly=(a,b) if s.shape[a]==Shape.DISC else (b,a)
                    vertices=world_polygon(s.polygons[poly],s.pos[poly],s.angle[poly])
                    gap=disc_polygon_separation(s.pos[disc],s.radius[disc],vertices)[0]
                if gap < deepest:
                    deepest, participants=gap,("pair",a,b)
            for wi, wall in enumerate(self.world.walls):
                if isinstance(wall, SegmentWall):
                    if self.wall_search == "bounded" and s.shape[a] != Shape.DISC:
                        nearest, _ = closest_point_segment(s.pos[a], wall.start, wall.end)
                        lower_bound = (float(np.linalg.norm(s.pos[a] - nearest)) -
                                       s.radius[a] - wall.thickness)
                        if lower_bound >= deepest + self.tol.geometry:
                            continue
                    if s.shape[a] == Shape.DISC:
                        nearest,_=closest_point_segment(s.pos[a],wall.start,wall.end)
                        gap=float(np.linalg.norm(s.pos[a]-nearest)-s.radius[a]-wall.thickness)
                    else:
                        vertices=world_polygon(s.polygons[a],s.pos[a],s.angle[a])
                        distance,*_=polygon_segment_witnesses(vertices,wall.start,wall.end)
                        gap=distance-wall.thickness
                    if gap < deepest:
                        deepest,participants=gap,("segment",a,wi)
                    continue
                if wall_points is not None:
                    lower_bound = (float(lower_bounds[a, wi]) if lower_bounds is not None
                                   else float((s.pos[a] - wall_points[wi]) @
                                              wall.inward_normal - s.radius[a]))
                    if lower_bound >= deepest + self.tol.geometry:
                        continue
                if lower_bounds is not None:
                    gap = (lower_bound if s.shape[a] == Shape.DISC else
                           polygon_wall_gap(s, a, wall_points[wi], wall.inward_normal))
                    if gap < deepest:
                        deepest, participants = gap, ("wall", a, wi)
                    continue
                if s.shape[a] == Shape.DISC:
                    support=s.pos[a]-s.radius[a]*wall.inward_normal
                else:
                    vertices=world_polygon(s.polygons[a],s.pos[a],s.angle[a])
                    support=vertices[np.argmin(vertices@wall.inward_normal)]
                point = wall_points[wi] if wall_points is not None else wall.point_at(self.time)
                gap=float((support-point)@wall.inward_normal)
                if gap < deepest:
                    deepest,participants=gap,("wall",a,wi)
        self.max_penetration=max(self.max_penetration,-deepest)
        return deepest,participants

    def _fail(self, reason: str, checkpoint: Checkpoint, **details: Any) -> None:
        self.failure_checkpoint=checkpoint
        self.failure_diagnostic={"reason":reason,"time":self.time,"seed":self.seed,
                                 "events":self.event_count,"details":details}
        raise NumericalFailure(f"{reason} at t={self.time}: {details}")

    def save_failure(self, directory: str | Path) -> Path:
        """Write the retained pre-failure state and diagnostic without pickle."""
        if self.failure_checkpoint is None or self.failure_diagnostic is None:
            raise ValueError("no numerical failure has been retained")
        destination=Path(directory)
        destination.mkdir(parents=True,exist_ok=True)
        cp=self.failure_checkpoint
        state=cp.bodies
        np.savez_compressed(destination/"state.npz",pos=state.pos,vel=state.vel,
                            angle=state.angle,omega=state.omega,mass=state.mass,
                            inertia=state.inertia,radius=state.radius,shape=state.shape,
                            species=state.species,dynamic=state.dynamic,version=state.version,
                            ids=state.ids,**{f"polygon_{i}":p for i,p in state.polygons.items()})
        mechanism=cp.mechanism
        shaft=getattr(mechanism,"shaft",None)
        payload={**self.failure_diagnostic,"rng_state":cp.rng_state,
                 "memberships":cp.memberships.tolist(),"event_count":cp.event_count,
                 "cluster_count":cp.cluster_count,
                 "max_cluster_residual":cp.max_cluster_residual,
                 "world_metadata":{k:v for k,v in self.world.metadata.items() if k!="clock"},
                 "mechanism":None if shaft is None else {
                     "phi":shaft.phi,"momentum":shaft.momentum,"inertia":shaft.inertia,
                     "prescribed_omega":shaft.prescribed_omega,"load_torque":shaft.load_torque}}
        (destination/"diagnostic.json").write_text(json.dumps(payload,indent=2,default=lambda x:
            x.tolist() if isinstance(x,np.ndarray) else x.item() if isinstance(x,np.generic) else str(x))+"\n")
        return destination

    def _advance_free(self, dt: float) -> None:
        old_positions = self.world.bodies.pos.copy() if self.world.portals else None
        self.world.bodies.advance(dt)
        if self.world.mechanism is not None:
            motor_before = getattr(self.world.mechanism,"motor_work",0.0)
            shaft=getattr(self.world.mechanism,"shaft",None)
            if self._controlled_phase_origin is not None and shaft is not None:
                target_phase=(self._controlled_phase_origin+
                              shaft.prescribed_omega*(self.time+dt))
                load_work=self.world.mechanism.advance(dt,target_phase=target_phase)
            else:
                load_work = self.world.mechanism.advance(dt)
            self.ledger.load_output.add(load_work)
            motor_work = getattr(self.world.mechanism,"motor_work",0.0)-motor_before
            self.ledger.work_on.add(motor_work-load_work)
        self.time += dt
        clock = self.world.metadata.get("clock")
        if clock is not None:
            clock[0] = self.time
        # Portal crossings are scheduled as physical events by _portal_toi.

    def _update_portals(self, old: np.ndarray) -> None:
        new = self.world.bodies.pos
        for portal in self.world.portals:
            n = portal.normal/np.linalg.norm(portal.normal)
            tangent = np.array([-n[1], n[0]])
            for i in range(self.world.bodies.n):
                a, b = float((old[i]-portal.point)@n), float((new[i]-portal.point)@n)
                if a*b > 0 or a == b:
                    continue
                u = a/(a-b)
                crossing = old[i]+u*(new[i]-old[i])
                if abs(float((crossing-portal.point)@tangent)) <= portal.half_length:
                    self.memberships[i] = portal.inside_label if b < 0 else 0
                    self.events.append(InteractionRecord(self.time-dt if False else self.time,
                        "portal_crossing", (i,), metadata={"host": portal.host_id,
                        "region": int(self.memberships[i])}))

    def _resolve(self, result: TOIResult) -> None:
        c = result.contact
        before = self.energy()
        gas_before = self.world.bodies.kinetic_energy()
        thermal_diagnostics=None
        if c.boundary is None:
            impulse, incoming = resolve_elastic(self.world.bodies, c)
            kind, heat, work = "elastic", 0.0, 0.0
            vector = tuple(impulse*np.asarray(c.normal))
            participants = (c.a, c.b)
        elif c.boundary is not None and c.boundary < 0:
            pi=-c.boundary-1; portal=self.world.portals[pi]
            direction=float(self.world.bodies.vel[c.a]@portal.normal)
            delta_u=portal.delta_u if direction < 0 else -portal.delta_u
            crossed, impulse=resolve_energy_step(self.world.bodies,c,delta_u)
            vector=tuple(impulse*np.asarray(c.normal)); participants=(c.a,)
            heat=work=0.0; kind="portal_crossing" if crossed else "portal_reflection"
            if crossed:
                # At the surface use post-event velocity to determine destination.
                direction=float(self.world.bodies.vel[c.a]@portal.normal)
                self.memberships[c.a]=portal.outside_label if direction > 0 else portal.inside_label
            self.events.append(InteractionRecord(self.time,kind,participants,vector,before,
                self.energy(),0.,0.,{"host":portal.host_id,"region":int(self.memberships[c.a]),
                "delta_u":delta_u}))
            self.event_count += 1
            return
        else:
            wall = self.world.walls[c.boundary]
            mechanism=getattr(wall,"mechanism",None)
            shaft=getattr(mechanism,"shaft",None)
            if wall.name=="cam_piston" and shaft is not None and shaft.prescribed_omega is None:
                n=np.asarray(c.normal); point=np.asarray(c.point); s=self.world.bodies
                particle_v=contact_velocity(s,c.a,point-s.pos[c.a])
                fp=mechanism.cam.piston_dx_dphi(shaft.phi)
                wall_v=np.array([fp*shaft.omega,0.])
                g=float((wall_v-particle_v)@n)
                d=inverse_effective_mass(s,c)+(fp*n[0])**2/shaft.effective_inertia()
                impulse=-2*g/d if g < 0 else 0.0
                apply_impulse(s,c,impulse)
                shaft.momentum += impulse*fp*n[0]
                heat=work=0.0; kind="cam_recoil"
            else:
                if wall.kind_at(self.time) in (BoundaryKind.HOT, BoundaryKind.COLD):
                    thermal_diagnostics={}
                impulse, heat, work = resolve_wall(self.world.bodies, c.a, wall,
                                                    np.asarray(c.point), self.rng, self.time,
                                                    np.asarray(c.normal),thermal_diagnostics)
                kind = wall.kind_at(self.time).value
                if (wall.name == "cam_piston" and shaft is not None and
                    shaft.prescribed_omega is not None):
                    mechanism.motor_work += work
            vector = tuple(impulse*np.asarray(c.normal))
            participants = (c.a,)
            if kind == "hot": self.ledger.heat_hot.add(heat)
            elif kind == "cold": self.ledger.heat_cold.add(heat)
            else: self.ledger.heat_other.add(heat)
            self.ledger.work_on.add(work)
            self.ledger.support_impulse_x.add(-vector[0])
            self.ledger.support_impulse_y.add(-vector[1])
        after = self.energy()
        if c.boundary is not None and c.boundary >= 0 and wall.name == "cam_piston":
            self.ledger.piston_work_on_gas.add(
                self.world.bodies.kinetic_energy()-gas_before)
        self.events.append(InteractionRecord(self.time, kind, participants, vector,
                           before, after, heat, work,
                           {"boundary": c.boundary,
                            **(thermal_diagnostics or {})}
                           if c.boundary is not None else None))
        self.event_count += 1

    def _resolve_contacts(self, result: TOIResult, checkpoint: Checkpoint) -> None:
        contacts = result.contacts or (result.contact,)
        if len(contacts) == 1:
            self._resolve(result)
            return
        groups: list[list[Contact]] = []
        for contact in contacts:
            matching = [group for group in groups if any(
                contact.a in (old.a, old.b) or
                (contact.b is not None and contact.b in (old.a, old.b))
                for old in group)]
            if not matching:
                groups.append([contact])
            else:
                matching[0].append(contact)
                for other in matching[1:]:
                    matching[0].extend(other)
                    groups.remove(other)
        for group in groups:
            if len(group) == 1:
                c=group[0]
                self._resolve(TOIResult(TOIStatus.COLLISION,result.time,result.error,c))
                continue
            for c in group:
                if c.boundary is None:
                    continue
                if c.boundary < 0:
                    self.restore(checkpoint)
                    self._fail("unsupported simultaneous portal contact",checkpoint,
                               contacts=[x.__dict__ for x in group])
                wall=self.world.walls[c.boundary]
                if (wall.kind_at(self.time) != BoundaryKind.SPECULAR or
                    np.linalg.norm(wall.velocity_at(self.time)) > self.tol.velocity or
                    (getattr(getattr(wall,"mechanism",None),"shaft",None) is not None and
                     wall.mechanism.shaft.prescribed_omega is None)):
                    self.restore(checkpoint)
                    self._fail("unsupported simultaneous moving or thermal contact",checkpoint,
                               contacts=[x.__dict__ for x in group])
            try:
                impulses,residual,rank=elastic_cluster_impulses(
                    self.world.bodies,tuple(group),self.tol.velocity)
            except ValueError as exc:
                self.restore(checkpoint)
                self._fail("ambiguous elastic cluster",checkpoint,
                           reason=str(exc),contacts=[x.__dict__ for x in group])
            before=self.energy()
            for c,impulse in zip(group,impulses):
                apply_impulse(self.world.bodies,c,float(impulse))
                if c.boundary is not None:
                    vector=impulse*np.asarray(c.normal)
                    self.ledger.support_impulse_x.add(-float(vector[0]))
                    self.ledger.support_impulse_y.add(-float(vector[1]))
            after=self.energy()
            if abs(after-before) > self.tol.energy*max(1.0,abs(before)):
                self.restore(checkpoint)
                self._fail("elastic cluster energy residual",checkpoint,
                           residual=after-before,contacts=[x.__dict__ for x in group])
            self.cluster_count += 1
            self.max_cluster_residual=max(self.max_cluster_residual,residual)
            for c,impulse in zip(group,impulses):
                vector=tuple(float(x) for x in impulse*np.asarray(c.normal))
                kind="elastic" if c.boundary is None else "specular"
                participants=(c.a,c.b) if c.b is not None else (c.a,)
                self.events.append(InteractionRecord(self.time,kind,participants,
                                   vector,before,after,0.0,0.0,
                                   {"boundary":c.boundary,"cluster_size":len(group),
                                    "cluster_rank":rank,"cluster_residual":residual}))
                self.event_count += 1

    def _prepare_interval(self) -> None:
        if self._pending is not None:
            return
        cp = self.checkpoint(include_history=False)
        horizon = self.max_horizon
        queries = []
        for _ in range(32):
            q = self._earliest(horizon)
            if q.status != TOIStatus.INDETERMINATE:
                break
            queries.append({"horizon": horizon, "reason": q.reason})
            self.ccd_refinements += 1
            horizon *= .5
            if horizon <= 16*self.tol.time:
                break
        if q.status == TOIStatus.INDETERMINATE:
            self.ccd_failures += 1
            self._fail("unresolved CCD", cp, queries=queries,
                       tolerance=self.tol.__dict__)
        mechanism = self.world.mechanism
        if mechanism is not None and hasattr(mechanism,"next_branch_transition"):
            transition = mechanism.next_branch_transition(horizon,self.tol.time)
            if transition is not None:
                branch_time, details = transition
                if self._controlled_phase_origin is not None:
                    sign=-1. if mechanism.reversed_cycle else 1.
                    absolute=(sign*details["phase"]-
                              self._controlled_phase_origin)/mechanism.shaft.prescribed_omega
                    branch_time=absolute-self.time
                if q.status == TOIStatus.NO_COLLISION or branch_time <= q.time+self.tol.time:
                    q = TOIResult(TOIStatus.BRANCH,branch_time,self.tol.time,
                                  transition=details)
        dt = horizon if q.status == TOIStatus.NO_COLLISION else float(q.time)
        self._pending = (cp, q, self.time+dt)

    def advance_to(self, physical_time: float) -> Snapshot:
        if physical_time < self.time-self.tol.time:
            raise ValueError("cannot run backward; restore a checkpoint instead")
        zero_events = 0
        while (self.time < physical_time-self.tol.time or
               (self._pending is not None and self._pending[2] <= physical_time+self.tol.time)):
            self._prepare_interval()
            cp, q, endpoint = self._pending
            # Materialize each observation directly from the same interval
            # origin. A caller's sampling cadence cannot perturb event search.
            self._restore_interval_origin(cp)
            if physical_time < endpoint-self.tol.time:
                self._advance_free(physical_time-cp.time)
                gap, participants = self._penetration()
                if gap < -self.tol.geometry:
                    self._restore_interval_origin(cp)
                    self._fail("penetration at sample",cp,gap=gap,participants=participants)
                self._pending = (cp, q, endpoint)
                break
            dt = endpoint-cp.time
            self._advance_free(dt)
            gap, participants=self._penetration()
            if gap < -self.tol.geometry:
                self.restore(cp)
                self._fail("penetration before impact",cp,gap=gap,participants=participants)
            if q.status == TOIStatus.COLLISION:
                self._resolve_contacts(q,cp)
                gap, participants=self._penetration()
                if gap < -self.tol.geometry:
                    self.restore(cp)
                    self._fail("penetration after impact",cp,gap=gap,participants=participants)
                zero_events = zero_events+1 if dt <= self.tol.time else 0
                if zero_events > 8:
                    self._fail("repeated zero-time contacts",cp)
            elif q.status == TOIStatus.BRANCH:
                mechanism=self.world.mechanism
                energy = self.energy()
                self.events.append(InteractionRecord(self.time,"branch_transition",(),
                                   energy_before=energy,energy_after=energy,
                                   metadata=q.transition))
                self.event_count += 1
                if q.transition["boundary"] == 0 and q.transition["to_branch"] == "hot":
                    direction=-1. if mechanism.reversed_cycle else 1.
                    if (not self.cycle_markers or direction*(mechanism.shaft.phi-
                        self.cycle_markers[-1].phase) >= 2*math.pi-1e-8):
                        self.cycle_markers.append(self._cycle_marker())
                zero_events = zero_events+1 if dt <= self.tol.time else 0
                if zero_events > 8:
                    self._fail("repeated zero-time branch transitions",cp)
            else:
                zero_events = 0
            self._pending = None
        snap = self.snapshot()
        return snap

    def step_collision(self) -> Snapshot:
        self._prepare_interval()
        while self._pending[1].status == TOIStatus.BRANCH:
            self.advance_to(self._pending[2])
            self._prepare_interval()
        return self.advance_to(self._pending[2])

    def snapshot(self) -> Snapshot:
        branch = self.world.mechanism.branch if self.world.mechanism is not None else None
        apparatus = (self.world.mechanism.apparatus_state()
                     if hasattr(self.world.mechanism,"apparatus_state") else ())
        phase = (self.world.mechanism.shaft.phi
                 if hasattr(self.world.mechanism,"shaft") else None)
        return Snapshot(self.time, self.world.bodies.pos.copy(), self.world.bodies.vel.copy(),
                        self.world.bodies.angle.copy(), self.world.bodies.omega.copy(),
                        self.energy(), translational_temperature(self.world.bodies),
                        rotational_temperature(self.world.bodies), branch,
                        tuple(np.bincount(self.memberships, minlength=4)), self.event_count,
                        self.ledger.residual(self.energy()), apparatus, phase)

    def apply_command(self, command: dict[str, Any]) -> dict[str, Any]:
        name = command.get("name")
        if name == "set_velocity":
            self._pending = None
            i = int(command["body"]); self.world.bodies.vel[i] = command["velocity"]
            self.world.bodies.version[i] += 1
            return {"ok": True, "time": self.time, "intervention": True}
        if name == "reverse_particle_velocities":
            self._pending = None
            moving=self.world.bodies.dynamic
            self.world.bodies.vel[moving] *= -1
            self.world.bodies.omega[moving] *= -1
            self.world.bodies.version[moving] += 1
            return {"ok": True, "time": self.time, "intervention": True}
        if name == "checkpoint":
            return {"ok": True, "checkpoint": self.checkpoint()}
        raise ValueError(f"unknown command {name!r}")
