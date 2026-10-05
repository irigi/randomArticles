from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from .base import Experiment
from .common import particle_specs
from ..config import RunConfig
from ..core.boundaries import BoundaryKind, Wall
from ..core.apparatus import ApparatusComponent
from ..core.mechanisms import CarnotCam, Shaft
from ..core.numeric import controlled_cam_disc_toi, controlled_cam_triangle_toi
from ..core.state import BodyArrays, Shape
from ..runner.simulation import World


@dataclass
class CarnotMechanism:
    cam: CarnotCam
    shaft: Shaft
    hot_temperature: float
    cold_temperature: float
    reversed_cycle: bool = False
    motor_work: float = 0.0
    load_output: float = 0.0

    @property
    def branch(self) -> str:
        return self.cam.branch(self.shaft.phi, self.reversed_cycle)

    def energy(self) -> float:
        return self.shaft.kinetic_energy()+self.shaft.spring_energy()

    def apparatus_state(self) -> tuple[ApparatusComponent, ...]:
        """Physical piston plus schematic, phase-linked drive components."""
        phi = self.shaft.phi
        height = self.cam.height
        piston_x = self.cam.piston_x(phi)
        piston_speed = self.cam.piston_dx_dphi(phi)*self.shaft.omega
        cylinder_width = max(self.cam.areas)/height
        shaft_x = .26*cylinder_width
        cy = -.48*height
        wheel_radius = .2*height
        cam_x = .53*cylinder_width
        selector_x = .8*cylinder_width
        sign = -1. if self.reversed_cycle else 1.

        def ring(cx: float, radius_at, rotation: float):
            points = []
            for index in range(129):
                local = 2*math.pi*index/128
                angle = local+rotation
                radius = radius_at(local)
                points.append((cx+radius*math.cos(angle),
                               cy+radius*math.sin(angle)))
            return tuple(points)

        min_x = min(self.cam.areas)/height
        span_x = (max(self.cam.areas)-min(self.cam.areas))/height
        cam_scale = .12*height/span_x if span_x else 0.
        piston_radius = lambda angle: (.15*height +
            cam_scale*(self.cam.piston_x(angle)-min_x))
        # World angle zero meets local phase phi because the track rotates -phi.
        follower_radius = piston_radius(phi)
        selector_levels = {"hot": .22, "adiabatic_expansion": .19,
                           "cold": .16, "adiabatic_compression": .19}
        selector_radius = lambda angle: (selector_levels[
            self.cam.branch(angle)]*height)
        branch = self.branch
        shoe_radius = selector_levels[branch]*height
        wheel = ring(shaft_x, lambda _: wheel_radius, phi)
        spring = tuple((shaft_x+(.06+.003*index)*height*
                        math.cos(phi-self.shaft.spring_rest+index*math.pi/4),
                        cy+(.06+.003*index)*height*
                        math.sin(phi-self.shaft.spring_rest+index*math.pi/4))
                       for index in range(25))
        piston_energy = .5*self.shaft.piston_mass*piston_speed**2
        flywheel_energy = .5*self.shaft.inertia*self.shaft.omega**2
        return (
            ApparatusComponent("cylinder", "outline",
                ((0.,0.),(piston_x,0.),(piston_x,height),(0.,height),(0.,0.))),
            ApparatusComponent("piston", "segment",
                ((piston_x,0.),(piston_x,height)),
                kinetic_energy=piston_energy,velocity=(piston_speed,0.)),
            ApparatusComponent("piston_cam", "closed_track",
                ring(cam_x,piston_radius,-phi),angular_velocity=-self.shaft.omega),
            ApparatusComponent("cam_follower", "point",
                ((cam_x+follower_radius,cy),),
                velocity=(self.cam.piston_dx_dphi(phi)*cam_scale*
                          self.shaft.omega,0.)),
            ApparatusComponent("selector_cam", "closed_track",
                ring(selector_x,selector_radius,-sign*phi),
                angular_velocity=-sign*self.shaft.omega,state=branch),
            ApparatusComponent("selector_shoe", "point",
                ((selector_x+shoe_radius,cy),),state=branch),
            ApparatusComponent("thermal_wall", "segment",
                ((0.,0.),(0.,height)),state=branch),
            ApparatusComponent("flywheel", "closed_track",wheel,
                kinetic_energy=flywheel_energy,angular_velocity=self.shaft.omega),
            ApparatusComponent("flywheel_spoke", "segment",
                ((shaft_x,cy),(shaft_x+wheel_radius*math.cos(phi),
                               cy+wheel_radius*math.sin(phi))),
                angular_velocity=self.shaft.omega),
            ApparatusComponent("shaft_spring", "polyline",spring,
                potential_energy=self.shaft.spring_energy(),
                state="torsion" if self.shaft.spring_k else "inactive"),
            ApparatusComponent("load", "segment",
                ((shaft_x+wheel_radius,cy),(shaft_x+wheel_radius+.25*height,cy)),
                work_output=self.load_output,angular_velocity=self.shaft.omega,
                state="opposing_torque" if self.shaft.load_torque else "inactive"),
        )

    def advance(self, dt: float, target_phase: float | None = None) -> float:
        before_energy = self.energy()
        if target_phase is not None and self.shaft.prescribed_omega is not None:
            travel=abs(target_phase-self.shaft.phi)
            self.shaft.phi=target_phase
        else:
            _,_,travel = self.shaft.advance_with_travel(dt)
        work = self.shaft.load_torque*travel
        if self.shaft.prescribed_omega is not None:
            self.motor_work += self.energy()-before_energy+work
        self.load_output += work
        return work

    def phi_at(self, absolute_time: float, now: float) -> float:
        return self.shaft.trajectory(max(0.0, absolute_time-now))[0]

    def omega_at(self, absolute_time: float, now: float) -> float:
        return self.shaft.trajectory(max(0.0, absolute_time-now))[2]

    def next_branch_transition(self, horizon: float, time_tol: float):
        """First cam sector crossing along the predicted shaft trajectory."""
        sign = -1.0 if self.reversed_cycle else 1.0
        period = 2*math.pi
        boundaries = self.cam.boundaries[:-1]
        knots = (0.0, *self.shaft.turning_times(horizon), horizon)
        for start, end in zip(knots[:-1],knots[1:]):
            a = sign*self.shaft.trajectory(start)[0]
            b = sign*self.shaft.trajectory(end)[0]
            if abs(b-a) <= self.shaft.speed_bound(end)*time_tol:
                continue
            direction = 1 if b > a else -1
            candidate = None
            for index, boundary in enumerate(boundaries):
                if direction > 0:
                    cycle = math.floor((a-boundary)/period)+1
                else:
                    cycle = math.ceil((a-boundary)/period)-1
                target = float(boundary+cycle*period)
                if abs(target-a) <= 4*self.shaft.speed_bound(horizon)*time_tol:
                    continue
                if (direction > 0 and target <= b or
                    direction < 0 and target >= b):
                    if candidate is None or direction*target < direction*candidate[0]:
                        candidate = (target,index)
            if candidate is None:
                continue
            target,index = candidate
            if self.shaft.prescribed_omega is not None:
                # Controlled phase is affine in time. Avoid accumulating a
                # one-sided bisection error at every quarter-cycle boundary.
                hi=(sign*target-self.shaft.phi)/self.shaft.prescribed_omega
            else:
                lo,hi = start,end
                for _ in range(80):
                    if hi-lo <= time_tol:
                        break
                    middle = .5*(lo+hi)
                    value = sign*self.shaft.trajectory(middle)[0]
                    if direction*value >= direction*target:
                        hi = middle
                    else:
                        lo = middle
            delta = 1e-8
            previous = self.cam.branch(target-direction*delta)
            following = self.cam.branch(target+direction*delta)
            return hi,{"phase":target,"boundary":index,
                       "from_branch":previous,"to_branch":following}
        return None


class CamPistonWall(Wall):
    def __init__(self, mechanism: CarnotMechanism, now_getter, height: float):
        super().__init__(np.array([mechanism.cam.piston_x(mechanism.shaft.phi), height/2]),
                         np.array([-1., 0.]), height, name="cam_piston")
        self.mechanism, self.now_getter = mechanism, now_getter

    def _phi(self, time):
        return self.mechanism.phi_at(time, self.now_getter())

    def point_at(self, time):
        return np.array([self.mechanism.cam.piston_x(self._phi(time)), self.point[1]])

    def velocity_at(self, time):
        phi = self._phi(time)
        return np.array([self.mechanism.cam.piston_dx_dphi(phi)*
                         self.mechanism.omega_at(time, self.now_getter()), 0.])

    def speed_bound(self, t0, t1):
        # Quintic derivative is bounded by 1.875 times the sector slope.
        b = self.mechanism.cam.boundaries
        slopes = [abs(self.mechanism.cam.areas[(i+1)%4]-self.mechanism.cam.areas[i])/(b[i+1]-b[i])
                  for i in range(4)]
        duration=max(0.0,t1-self.now_getter())
        return (1.875*max(slopes)*self.mechanism.shaft.speed_bound(duration)
                /self.mechanism.cam.height)

    def batch_disc_toi(self, state, possible, wall_speed, horizon, tol):
        if self.mechanism.shaft.prescribed_omega is None:
            return None
        return controlled_cam_disc_toi(state, possible, self.mechanism,
                                       wall_speed, horizon, tol)

    def batch_triangle_toi(self, state, possible, wall_speed, horizon, tol):
        if self.mechanism.shaft.prescribed_omega is None:
            return None
        return controlled_cam_triangle_toi(state, possible, self.mechanism,
                                           wall_speed, horizon, tol)


class SelectorWall(Wall):
    def __init__(self, mechanism: CarnotMechanism, height: float, now_getter):
        super().__init__(np.array([0., height/2]), np.array([1., 0.]), height, name="selector")
        self.mechanism, self.now_getter = mechanism, now_getter

    def kind_at(self, time):
        phi = self.mechanism.phi_at(time, self.now_getter())
        branch = self.mechanism.cam.branch(phi, self.mechanism.reversed_cycle)
        if branch == "hot": return BoundaryKind.HOT
        if branch == "cold": return BoundaryKind.COLD
        return BoundaryKind.SPECULAR

    def temperature_at(self, time):
        kind = self.kind_at(time)
        if kind == BoundaryKind.HOT: return self.mechanism.hot_temperature
        if kind == BoundaryKind.COLD: return self.mechanism.cold_temperature
        return None


class ThermalJacketWall(Wall):
    """Stationary exchanger, thermal in enabled sectors and specular otherwise."""

    def __init__(self, point, normal, length, name, mechanism, now_getter,
                 hot: bool, cold: bool):
        super().__init__(np.asarray(point), np.asarray(normal), length, name=name)
        self.mechanism = mechanism
        self.now_getter = now_getter
        self.hot = hot
        self.cold = cold

    def kind_at(self, time):
        phi = self.mechanism.phi_at(time, self.now_getter())
        branch = self.mechanism.cam.branch(phi, self.mechanism.reversed_cycle)
        if branch == "hot" and self.hot:
            return BoundaryKind.HOT
        if branch == "cold" and self.cold:
            return BoundaryKind.COLD
        return BoundaryKind.SPECULAR

    def temperature_at(self, time):
        kind = self.kind_at(time)
        if kind == BoundaryKind.HOT:
            return self.mechanism.hot_temperature
        if kind == BoundaryKind.COLD:
            return self.mechanism.cold_temperature
        return None


class ColdJacketWall(ThermalJacketWall):
    """Compatibility wrapper for the original cold-only jacket."""

    def __init__(self, point, normal, length, name, mechanism, now_getter):
        super().__init__(point, normal, length, name, mechanism, now_getter,
                         hot=False, cold=True)


class CarnotExperiment(Experiment):
    def __init__(self, triangles: bool = False):
        self.triangles = triangles

    def build(self, config: RunConfig, rng: np.random.Generator) -> World:
        height, a1, ratio = 1.0, 1.25, 1.35
        th = 1.5*config.temperature
        tc = th/config.temperature_ratio
        dof = 3 if self.triangles else 2
        base_cam = CarnotCam.design(a1, ratio, th, tc, dof, height)
        cam = CarnotCam(base_cam.areas, tuple(config.cam_fractions), height)
        omega = (-1 if config.reversed_cycle else 1)*config.shaft_speed
        prescribed = omega if config.shaft_mode == "controlled" else None
        shaft = Shaft(0.0, 20.0*omega, 20.0, prescribed_omega=prescribed,
                      load_torque=0.02, piston_mass=0.5, cam=cam)
        mech = CarnotMechanism(cam, shaft, th, tc, config.reversed_cycle)
        shape = Shape.TRIANGLE if self.triangles else Shape.DISC
        # Keep the occupied area fraction fixed as particle count changes.
        # At 48 particles triangles fill about 3.1% of the minimum cylinder.
        radius = 0.025 * math.sqrt(48 / config.particles) * config.carnot_radius_scale
        specs = particle_specs(config.particles, max(cam.areas)/height, height, radius,
                               (config.temperature if config.initial_temperature is None
                                else config.initial_temperature),
                               rng, shape, x_max=min(cam.areas)/height)
        state = BodyArrays.from_specs(specs)
        # The simulation updates this clock at every committed or sampled time.
        clock = [0.0]
        bottom_point, bottom_normal = np.array([1., 0.]), np.array([0., 1.])
        top_point, top_normal = np.array([1., height]), np.array([0., -1.])
        if config.cold_jacket and not config.hot_jacket:
            bottom = ColdJacketWall(bottom_point, bottom_normal, 3., "bottom",
                                    mech, lambda: clock[0])
            top = ColdJacketWall(top_point, top_normal, 3., "top",
                                 mech, lambda: clock[0])
        elif config.cold_jacket or config.hot_jacket:
            bottom = ThermalJacketWall(bottom_point, bottom_normal, 3., "bottom",
                                       mech, lambda: clock[0], config.hot_jacket,
                                       config.cold_jacket)
            top = ThermalJacketWall(top_point, top_normal, 3., "top",
                                    mech, lambda: clock[0], config.hot_jacket,
                                    config.cold_jacket)
        else:
            bottom = Wall(bottom_point, bottom_normal, 3., name="bottom")
            top = Wall(top_point, top_normal, 3., name="top")
        walls: list[Wall] = [SelectorWall(mech, height, lambda: clock[0]),
                             bottom, top, CamPistonWall(mech, lambda: clock[0], height)]
        return World(state, walls, mechanism=mech,
                     metadata={"name": "carnot_triangles" if self.triangles else "carnot_discs",
                               "clock": clock, "T_hot": th, "T_cold": tc,
                               "shaft_mode": config.shaft_mode,
                               "piston_mass": shaft.piston_mass,
                               "cold_jacket": config.cold_jacket,
                               "hot_jacket": config.hot_jacket,
                               "cam_fractions": cam.fractions,
                               "particle_radius": radius,
                               "carnot_radius_scale": config.carnot_radius_scale,
                               "particle_area_fraction_min": (config.particles * radius**2 *
                                   (3*math.sqrt(3)/4 if self.triangles else math.pi) / a1),
                               "ideal_efficiency": 1-tc/th, "cam_areas": cam.areas})
