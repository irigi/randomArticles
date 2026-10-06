"""`Simulation`-compatible driver for the compiled event-driven disc kernel.

The worker, instruments, exports and replay writer call the same methods and
attributes as on the reference `Simulation`. Supported worlds: smooth
spinless discs inside an axis-aligned box of stationary infinite-line walls
(specular or thermal), optionally with fixed circular posts, and with no
portals or mechanism.
"""
from __future__ import annotations

from dataclasses import dataclass
import copy
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from ..core import edmd
from ..core.boundaries import BoundaryKind, CircularPost, Wall
from ..core.events import InteractionRecord
from ..core.state import Shape, Tolerances
from ..measurements.ledger import CompensatedCounter, EnergyLedger
from ..measurements.observables import rotational_temperature, translational_temperature
from .simulation import NumericalFailure, Snapshot, World

_LOG_KINDS = ("elastic", "specular", "hot", "cold", "post")
_WALL_KINDS = {BoundaryKind.SPECULAR: edmd.WALL_SPECULAR,
               BoundaryKind.HOT: edmd.WALL_HOT, BoundaryKind.COLD: edmd.WALL_COLD}
_ARRAYS = ("pos", "vel", "tl", "cnt", "cell", "head", "nxt", "prv",
           "ht", "hk", "si", "sf", "rand", "pimp")
_GAP_KINDS = {edmd.EV_PAIR: "pair", edmd.EV_WALL: "wall", edmd.EV_POST: "post",
              edmd.EV_CELL: "cell"}


@dataclass
class EdmdCheckpoint:
    time: float
    arrays: dict[str, np.ndarray]
    rng_state: dict
    initial_energy: float
    history: tuple[InteractionRecord, ...] | None
    history_start: int


def edmd_unsupported_reason(world: World) -> str | None:
    """None when the world fits the kernel, otherwise why it does not."""
    s = world.bodies
    if world.mechanism is not None:
        return "mechanisms are not supported"
    if world.portals:
        return "portals are not supported"
    if np.any(s.shape != Shape.DISC) or np.any(np.isfinite(s.inertia)):
        return "only smooth spinless discs are supported"
    if not np.all(s.dynamic):
        return "static bodies are not supported"
    for wall in world.walls:
        if type(wall) is not Wall:
            return f"wall {wall.name!r} of type {type(wall).__name__} is not supported"
        if np.any(wall.velocity != 0):
            return f"wall {wall.name!r} moves"
    box = _box_bounds(world.walls)
    if box is None:
        return "walls must form an axis-aligned box"
    for post in world.posts:
        if not isinstance(post, CircularPost):
            return f"unsupported post {post!r}"
        (x, y), r = post.center, post.radius
        if x + r <= box[0] or x - r >= box[2] or y + r <= box[1] or y - r >= box[3]:
            return f"post {post.name!r} lies outside the wall box"
    return None


def _box_bounds(walls: list[Wall]) -> tuple[float, float, float, float] | None:
    bounds: dict[tuple[int, int], float] = {}
    for wall in walls:
        n = wall.inward_normal
        key = (int(round(n[0])), int(round(n[1])))
        if abs(abs(n[0]) + abs(n[1]) - 1) > 1e-12 or key in bounds:
            return None
        bounds[key] = float(wall.point[0] if key[1] == 0 else wall.point[1])
    try:
        box = bounds[(1, 0)], bounds[(0, 1)], bounds[(-1, 0)], bounds[(0, -1)]
    except KeyError:
        return None
    return box if box[2] > box[0] and box[3] > box[1] else None


class EdmdSimulation:
    """Compiled event-driven engine with the reference engine's interface."""

    numeric_backend = "edmd"

    def __init__(self, world: World, seed: int = 123,
                 tolerances: Tolerances = Tolerances(), record_events: bool = True,
                 random_block: int = 4096, log_capacity: int = 1 << 16):
        if not edmd.numba_available():
            raise RuntimeError("the edmd engine requires microthermo[accel] (numba)")
        reason = edmd_unsupported_reason(world)
        if reason is not None:
            raise ValueError(f"edmd engine: {reason}")
        self.world, self.seed, self.tol = world, seed, tolerances
        self.rng = np.random.default_rng(seed)
        self.random_block = random_block
        s = world.bodies
        n = s.n
        self.pos, self.vel = s.pos.copy(), s.vel.copy()
        self.tl = np.zeros(n)
        self.rad, self.mass = s.radius.copy(), s.mass.copy()
        self.cnt = np.zeros(n, np.int64)
        walls = world.walls
        self.wpt = np.asarray([w.point for w in walls], float).reshape(-1, 2)
        self.wnrm = np.asarray([w.inward_normal for w in walls], float).reshape(-1, 2)
        self.wkind = np.asarray([_WALL_KINDS[w.kind] for w in walls], np.int64)
        self.wtemp = np.asarray([w.temperature or 0.0 for w in walls], float)
        posts = world.posts
        self.ppos = np.asarray([p.center for p in posts], float).reshape(-1, 2)
        self.prad = np.asarray([p.radius for p in posts], float)
        self.pimp = np.zeros(4*len(posts))
        x0, y0, x1, y1 = _box_bounds(walls)
        side = 2.0*float(np.max(self.rad))
        nx = max(1, int((x1-x0)/side))
        ny = max(1, int((y1-y0)/side))
        self.fp = np.array([tolerances.geometry, tolerances.velocity, x0, y0,
                            (x1-x0)/nx, (y1-y0)/ny])
        self.ip = np.array([nx, ny], np.int64)
        self.cps, self.cpi = edmd.post_cell_lists(
            self.ppos, self.prad, float(np.max(self.rad)) + 1e-6*max(x1-x0, y1-y0),
            self.fp, self.ip)
        self.cell = np.zeros(n, np.int64)
        self.head = np.full(nx*ny, -1, np.int64)
        self.nxt = np.full(n, -1, np.int64)
        self.prv = np.full(n, -1, np.int64)
        capacity = 32*n + 4*(n + len(walls) + len(posts) + 2) + 1024
        self.ht = np.zeros(capacity)
        self.hk = np.zeros((capacity, 5), np.int64)
        self.si = np.zeros(edmd.SI_SIZE, np.int64)
        self.si[edmd.SI_LOG_ON] = 1 if record_events else 0
        self.sf = np.zeros(edmd.SF_SIZE)
        self.rand = np.zeros(0)
        self.lt = np.zeros(log_capacity)
        self.lk = np.zeros((log_capacity, 4), np.int64)
        self.lf = np.zeros((log_capacity, 7))
        self._log_chunks: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
        self._events: list[InteractionRecord] = []
        self._event_history_start = 0
        self.samples: list[Snapshot] = []
        self.memberships = np.zeros(n, dtype=np.int16)
        self.max_penetration = 0.0
        self.ccd_refinements = self.ccd_failures = self.cluster_count = 0
        self.max_cluster_residual = 0.0
        self.cycle_markers: list = []
        self.failure_checkpoint: EdmdCheckpoint | None = None
        self.failure_diagnostic: dict[str, Any] | None = None
        bad = edmd.build_grid(self._bodies(), self._grid(), 0.0)
        if bad >= 0:
            raise ValueError(f"body {bad} lies outside the wall box")
        gap, *_ = self._min_gap(0.0)
        if gap < -tolerances.geometry:
            raise ValueError(f"invalid overlapping initial state, gap={gap}")
        self.sf[edmd.SF_KE] = s.kinetic_energy()
        self._initial_energy = s.kinetic_energy()
        self._rebuild()

    # -- kernel plumbing --------------------------------------------------

    # Tuples are rebuilt on every call because restore() and heap growth
    # replace the arrays.
    def _bodies(self):
        return (self.pos, self.vel, self.tl, self.rad, self.mass, self.cnt)

    def _grid(self):
        return (self.cell, self.head, self.nxt, self.prv, self.fp, self.ip)

    def _static(self):
        return (self.wpt, self.wnrm, self.wkind, self.wtemp,
                self.ppos, self.prad, self.cps, self.cpi)

    def _calendar(self):
        return (self.ht, self.hk, self.si, self.sf)

    def _io(self):
        return (self.rand, self.lt, self.lk, self.lf, self.pimp)

    def _rebuild(self) -> None:
        if not edmd.rebuild(self._bodies(), self._grid(), self._static(), self._calendar()):
            self._overlap_failure("overlap while predicting events")

    def _min_gap(self, t: float):
        return edmd.min_gap(t, self._bodies(), self._grid(), self._static())

    def _run(self, t_end: float, max_events: int) -> int:
        """Run the kernel, servicing buffer, log and heap requests."""
        processed = 0
        while True:
            before = int(self.si[edmd.SI_EVENTS])
            status = edmd.advance(t_end, max_events - processed, self._bodies(),
                                  self._grid(), self._static(), self._calendar(),
                                  self._io())
            processed += int(self.si[edmd.SI_EVENTS]) - before
            if status == edmd.NEED_RANDOM:
                self.rand = self.rng.random(self.random_block)
                self.si[edmd.SI_RAND] = 0
            elif status == edmd.LOG_FULL:
                self._flush_log()
            elif status == edmd.HEAP_FULL:
                capacity = 2*len(self.ht)
                self.ht = np.zeros(capacity)
                self.hk = np.zeros((capacity, 5), np.int64)
                self._rebuild()
            elif status in (edmd.ERR_OVERLAP, edmd.ERR_CELL):
                self._flush_log()
                self._overlap_failure("overlap during event prediction"
                                      if status == edmd.ERR_OVERLAP else
                                      "body left the cell grid")
            else:
                self._flush_log()
                return processed

    def _flush_log(self) -> None:
        k = int(self.si[edmd.SI_LOG])
        if k:
            self._log_chunks.append((self.lt[:k].copy(), self.lk[:k].copy(),
                                     self.lf[:k].copy()))
            self.si[edmd.SI_LOG] = 0

    def _materialize_events(self) -> None:
        for lt, lk, lf in self._log_chunks:
            for t, (kind, a, b, wall), f in zip(lt.tolist(), lk.tolist(), lf.tolist()):
                name = _LOG_KINDS[kind]
                if kind == edmd.LOG_ELASTIC:
                    self._events.append(InteractionRecord(
                        t, name, (a, b), (f[0], f[1]), f[3], f[4]))
                    continue
                if kind == edmd.LOG_POST:
                    self._events.append(InteractionRecord(
                        t, name, (a,), (f[0], f[1]), f[3], f[4], 0.0, 0.0,
                        {"post": wall}))
                    continue
                metadata: dict[str, Any] = {"boundary": wall}
                if kind != edmd.LOG_SPECULAR:
                    metadata.update(reservoir_temperature=float(self.wtemp[wall]),
                                    incoming_normal_mode_energy=f[5],
                                    outgoing_normal_mode_energy=f[6])
                self._events.append(InteractionRecord(
                    t, name, (a,), (f[0], f[1]), f[3], f[4], f[2], 0.0, metadata))
        self._log_chunks.clear()

    def _sync_world(self) -> None:
        edmd.positions_at(self.time, self.pos, self.vel, self.tl, self.world.bodies.pos)
        self.world.bodies.vel[:] = self.vel

    def _check_gap(self, reason: str) -> None:
        gap, kind, a, b = self._min_gap(self.time)
        if gap < 0:
            self.max_penetration = max(self.max_penetration, -gap)
        if gap < -self.tol.geometry:
            self._fail(reason, gap=gap, participants=(_GAP_KINDS[kind], int(a), int(b)))

    def _overlap_failure(self, reason: str) -> None:
        si = self.si
        self._fail(reason, gap=float(self.sf[edmd.SF_ERR_GAP]),
                   participants=(_GAP_KINDS[int(si[edmd.SI_ERR_KIND])],
                                 int(si[edmd.SI_ERR_A]), int(si[edmd.SI_ERR_B])))

    def _fail(self, reason: str, **details: Any) -> None:
        self.failure_checkpoint = self.checkpoint(include_history=False)
        self.failure_diagnostic = {"reason": reason, "time": self.time, "seed": self.seed,
                                   "events": self.event_count, "details": details}
        raise NumericalFailure(f"{reason} at t={self.time}: {details}")

    # -- public interface (mirrors runner.simulation.Simulation) ----------

    @property
    def time(self) -> float:
        return float(self.sf[edmd.SF_TIME])

    @property
    def event_count(self) -> int:
        return int(self.si[edmd.SI_EVENTS])

    @property
    def events(self) -> list[InteractionRecord]:
        self._materialize_events()
        return self._events

    @events.setter
    def events(self, value: list[InteractionRecord]) -> None:
        self._log_chunks.clear()
        self._events = value

    def drain_events(self) -> list[InteractionRecord]:
        drained = self.events
        self._events = []
        self._event_history_start = self.event_count
        return drained

    @property
    def ledger(self) -> EnergyLedger:
        sf = self.sf

        def counter(k: int) -> CompensatedCounter:
            return CompensatedCounter(float(sf[k]), float(sf[k+1]))

        return EnergyLedger(self._initial_energy, counter(edmd.SF_HEAT_HOT),
                            counter(edmd.SF_HEAT_COLD), counter(edmd.SF_HEAT_OTHER),
                            support_impulse_x=counter(edmd.SF_SUPPORT_X),
                            support_impulse_y=counter(edmd.SF_SUPPORT_Y))

    def post_impulses(self) -> np.ndarray:
        """Total impulse delivered to each post so far, shape (posts, 2)."""
        return self.pimp.reshape(-1, 4)[:, [0, 2]].copy()

    def energy(self) -> float:
        return float(0.5*np.sum(self.mass*np.sum(self.vel*self.vel, axis=1)))

    def advance_to(self, physical_time: float) -> Snapshot:
        if physical_time < self.time - self.tol.time:
            raise ValueError("cannot run backward; restore a checkpoint instead")
        self._run(float(physical_time), np.iinfo(np.int64).max)
        self.sf[edmd.SF_TIME] = max(self.time, float(physical_time))
        self._check_gap("penetration at sample")
        self._sync_world()
        return self.snapshot()

    def step_collision(self) -> Snapshot:
        self._run(math.inf, 1)
        self._sync_world()
        return self.snapshot()

    def snapshot(self) -> Snapshot:
        s = self.world.bodies
        energy = self.energy()
        return Snapshot(self.time, s.pos.copy(), s.vel.copy(), s.angle.copy(),
                        s.omega.copy(), energy, translational_temperature(s),
                        rotational_temperature(s), None,
                        tuple(np.bincount(self.memberships, minlength=4)),
                        self.event_count, self.ledger.residual(energy))

    def checkpoint(self, include_history: bool = True) -> EdmdCheckpoint:
        return EdmdCheckpoint(
            self.time, {name: getattr(self, name).copy() for name in _ARRAYS},
            copy.deepcopy(self.rng.bit_generator.state), self._initial_energy,
            tuple(copy.deepcopy(self.events)) if include_history else None,
            self._event_history_start)

    def restore(self, cp: EdmdCheckpoint) -> None:
        for name, value in cp.arrays.items():
            setattr(self, name, value.copy())
        self.rng.bit_generator.state = copy.deepcopy(cp.rng_state)
        self._initial_energy = cp.initial_energy
        self._log_chunks.clear()
        if cp.history is not None:
            self._events = list(copy.deepcopy(cp.history))
            self._event_history_start = self.event_count - len(self._events)
        else:
            retained = self.event_count - self._event_history_start
            if retained < 0:
                raise ValueError("cannot restore discarded event history")
            self._events = self._events[:retained]
        self._sync_world()

    def save_failure(self, directory: str | Path) -> Path:
        if self.failure_checkpoint is None or self.failure_diagnostic is None:
            raise ValueError("no numerical failure has been retained")
        destination = Path(directory)
        destination.mkdir(parents=True, exist_ok=True)
        cp = self.failure_checkpoint
        np.savez_compressed(destination/"state.npz", **cp.arrays, radius=self.rad,
                            mass=self.mass, wall_points=self.wpt, wall_normals=self.wnrm,
                            wall_kinds=self.wkind, wall_temperatures=self.wtemp,
                            post_centers=self.ppos, post_radii=self.prad)
        payload = {**self.failure_diagnostic, "engine": "edmd", "rng_state": cp.rng_state,
                   "world_metadata": self.world.metadata}
        (destination/"diagnostic.json").write_text(json.dumps(
            payload, indent=2, default=lambda x: x.tolist() if isinstance(x, np.ndarray)
            else x.item() if isinstance(x, np.generic) else str(x)) + "\n")
        return destination

    def _bring_to_now(self) -> None:
        """Make every stored position current; call before changing velocities."""
        edmd.positions_at(self.time, self.pos, self.vel, self.tl, self.pos)
        self.tl[:] = self.time

    def _changed_velocities(self) -> None:
        """Replan after an intervention changed velocities at the current time."""
        self.cnt += 1
        self.sf[edmd.SF_KE] = self.energy()
        self._rebuild()
        self._sync_world()

    def apply_command(self, command: dict[str, Any]) -> dict[str, Any]:
        name = command.get("name")
        if name == "set_velocity":
            i = int(command["body"])
            self._bring_to_now()
            self.vel[i] = command["velocity"]
            self._changed_velocities()
            return {"ok": True, "time": self.time, "intervention": True}
        if name == "reverse_particle_velocities":
            self._bring_to_now()
            self.vel *= -1
            self._changed_velocities()
            return {"ok": True, "time": self.time, "intervention": True}
        if name == "checkpoint":
            return {"ok": True, "checkpoint": self.checkpoint()}
        raise ValueError(f"unknown command {name!r}")
