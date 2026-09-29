from __future__ import annotations

from dataclasses import dataclass, asdict
import math
import numpy as np

from .core.boundaries import BoundaryKind, Wall, resolve_wall
from .core.ccd import body_pair_toi, disc_disc_toi, disc_wall_toi
from .core.contacts import relative_normal_velocity, resolve_elastic, resolve_energy_step
from .core.events import Contact, TOIStatus
from .core.geometry import aperture_clearance, convex_separation, world_polygon
from .core.mechanisms import CarnotCam, LinearSpring
from .core.state import BodyArrays, BodySpec, Shape, Tolerances


@dataclass
class ValidationResult:
    name: str
    passed: bool
    metric: float
    tolerance: float
    detail: str = ""


def _result(name, metric, tolerance, detail=""):
    return ValidationResult(name, bool(metric <= tolerance), float(metric), float(tolerance), detail)


def run_validation(suite: str = "scientific") -> list[ValidationResult]:
    tol = Tolerances()
    results: list[ValidationResult] = []
    q = disc_disc_toi(np.array([0.,0.]), np.array([1.,0.]), .5,
                      np.array([3.,0.]), np.array([-1.,0.]), .5, 2., tol)
    err = abs((q.time if q.time is not None else 99)-1.0)
    results.append(_result("disc_disc_toi", err, 1e-12, q.status.value))
    tangent = disc_disc_toi(np.array([0.,0.]), np.array([1.,0.]), .5,
                            np.array([1.,1.]), np.array([0.,0.]), .5, 2., tol)
    results.append(_result("disc_tangency", abs((tangent.time or 0)-1.0), 1e-12, tangent.status.value))
    wall = disc_wall_toi(np.array([1.,0.]), np.array([-2.,0.]), .2,
                         np.array([0.,0.]), np.array([1.,0.]), 1.)
    results.append(_result("disc_wall_toi", abs((wall.time or 0)-.4), 1e-12, wall.status.value))

    rng = np.random.default_rng(77)
    worst_e = worst_p = worst_l = 0.0
    for _ in range(500):
        specs = [BodySpec(tuple(rng.normal(size=2)), tuple(rng.normal(size=2)),
                          mass=float(rng.uniform(.2,3)), radius=.1),
                 BodySpec(tuple(rng.normal(size=2)), tuple(rng.normal(size=2)),
                          mass=float(rng.uniform(.2,3)), radius=.1)]
        state = BodyArrays.from_specs(specs)
        delta = state.pos[1]-state.pos[0]
        n = delta/np.linalg.norm(delta)
        # Force an approaching normal relative velocity.
        if (state.vel[1]-state.vel[0])@n >= 0: state.vel[1] -= 2*n
        point = .5*(state.pos[0]+state.pos[1])
        contact = Contact(0, 1, tuple(point), tuple(n))
        e0, p0, l0 = state.kinetic_energy(), state.momentum(), state.angular_momentum()
        resolve_elastic(state, contact)
        scale = max(1., abs(e0))
        worst_e = max(worst_e, abs(state.kinetic_energy()-e0)/scale)
        worst_p = max(worst_p, np.linalg.norm(state.momentum()-p0)/max(1.,np.linalg.norm(p0)))
        worst_l = max(worst_l, abs(state.angular_momentum()-l0)/max(1.,abs(l0)))
    results += [_result("binary_energy", worst_e, 2e-14),
                _result("binary_momentum", worst_p, 2e-14),
                _result("binary_angular_momentum", worst_l, 2e-14)]

    # Exercise real polygon CCD witnesses and the angular impulse together.
    triangle_rng=np.random.default_rng(2026)
    collisions=0
    worst_e=worst_p=worst_l=worst_gap=0.0
    for _ in range(200):
        specs=[BodySpec((-.5,float(triangle_rng.uniform(-.15,.15))),(1.,0.),
                        mass=float(triangle_rng.uniform(.5,2.)),radius=.2,
                        shape=Shape.TRIANGLE,
                        angle=float(triangle_rng.uniform(-math.pi,math.pi)),
                        omega=float(triangle_rng.uniform(-1.,1.))),
               BodySpec((.5,float(triangle_rng.uniform(-.15,.15))),(-1.,0.),
                        mass=float(triangle_rng.uniform(.5,2.)),radius=.2,
                        shape=Shape.TRIANGLE,
                        angle=float(triangle_rng.uniform(-math.pi,math.pi)),
                        omega=float(triangle_rng.uniform(-1.,1.)))]
        state=BodyArrays.from_specs(specs)
        query=body_pair_toi(state,0,1,1.,tol)
        if query.status != TOIStatus.COLLISION:
            continue
        collisions+=1
        state.advance(query.time)
        if relative_normal_velocity(state,query.contact) >= 0.0:
            worst_gap=math.inf
            continue
        e0,p0,l0=state.kinetic_energy(),state.momentum(),state.angular_momentum()
        resolve_elastic(state,query.contact)
        worst_e=max(worst_e,abs(state.kinetic_energy()-e0)/max(1.,e0))
        worst_p=max(worst_p,float(np.linalg.norm(state.momentum()-p0))/max(1.,float(np.linalg.norm(p0))))
        worst_l=max(worst_l,abs(state.angular_momentum()-l0)/max(1.,abs(l0)))
        pa=world_polygon(state.polygons[0],state.pos[0],state.angle[0])
        pb=world_polygon(state.polygons[1],state.pos[1],state.angle[1])
        worst_gap=max(worst_gap,max(0.,-convex_separation(pa,pb)[0]))
    results += [_result("triangle_impact_coverage",200-collisions,10,
                        f"{collisions} of 200 approaching scenes collided"),
                _result("triangle_binary_energy",worst_e,2e-14),
                _result("triangle_binary_momentum",worst_p,2e-14),
                _result("triangle_binary_angular_momentum",worst_l,2e-14),
                _result("triangle_impact_penetration",worst_gap,tol.geometry)]

    state = BodyArrays.from_specs([BodySpec((0.,0.), (2.,0.), mass=1., radius=.1),
                                   BodySpec((1.,0.), (0.,0.), mass=4., radius=.1)])
    c = Contact(0, 1, (.5,0.), (1.,0.))
    total0 = state.kinetic_energy()
    crossed, _ = resolve_energy_step(state, c, .3)
    total1 = state.kinetic_energy()+(.3 if crossed else 0.)
    results.append(_result("portal_energy_step", abs(total1-total0), 2e-14, f"crossed={crossed}"))

    osc = LinearSpring(.7, -.2, 2., 3., .1)
    e0 = osc.energy()
    for _ in range(10000): osc.advance_exact(.01)
    results.append(_result("spring_exact_flow", abs(osc.energy()-e0)/e0, 2e-12))

    cam = CarnotCam.design(1., 1.4, 2., 1., 3)
    continuity = max(abs(cam.area(2*np.pi*i/4-1e-8)-cam.area(2*np.pi*i/4+1e-8)) for i in range(4))
    results.append(_result("cam_continuity", continuity, 1e-10))
    results.append(ValidationResult("aperture_clearance", aperture_clearance(.21,.1) > 0 and aperture_clearance(.19,.1) < 0,
                                    0., 0., "wide passes; narrow blocks"))

    # Thermal normal mode: transformed variable y=m v_n^2/(2T) has mean one.
    thermal = BodyArrays.from_specs([BodySpec((.1,0.), (-1.,0.), radius=.1)])
    tw = Wall(np.array([0.,0.]), np.array([1.,0.]), 1., BoundaryKind.HOT, 1.7)
    samples=[]
    for _ in range(20000):
        thermal.vel[0]=[-1.,0.]
        resolve_wall(thermal, 0, tw, np.array([0.,0.]), rng)
        samples.append(thermal.mass[0]*thermal.vel[0,0]**2/(2*1.7))
    results.append(_result("thermal_flux_mean", abs(float(np.mean(samples))-1.), .025))
    return results


def results_as_dict(results):
    return [asdict(r) for r in results]
