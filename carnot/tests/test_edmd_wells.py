"""Binding wells inside ring hosts (osmosis plan, milestone 4)."""

import math
import unittest

import numpy as np

from microthermo.core import edmd
from microthermo.core.boundaries import BoundaryKind
from microthermo.core.state import BodyArrays, BodySpec, RingGeometry, Shape
from microthermo.experiments.common import box_walls
from microthermo.runner.edmd_simulation import EdmdSimulation, edmd_unsupported_reason
from microthermo.runner.simulation import World

WELL = RingGeometry(.12, .15, .08, well_radius=.08, well_depth=1.5)


def disc(x, y, vx=0., vy=0., omega=0., r=.02):
    return BodySpec((x, y), (vx, vy), radius=r, omega=omega, inertia=.5*r*r)


def host(x, y, vx=0., vy=0., angle=0., omega=0., m=25., ring=WELL, inertia=None):
    return BodySpec((x, y), (vx, vy), mass=m, radius=ring.outer_radius, angle=angle,
                    omega=omega, inertia=ring.uniform_inertia(m) if inertia is None else inertia,
                    shape=Shape.HOST)


def make(specs, rings, width=2., height=1., temperature=None, roughness=1.):
    walls = box_walls(width, height)
    if temperature is not None:
        for w in walls:
            w.kind, w.temperature = BoundaryKind.HOT, temperature
    return EdmdSimulation(World(BodyArrays.from_specs(specs), walls, rings=rings,
                                contact_roughness=roughness))


def occupancy_box(depth, temperature, discs=60, seed=1, radius=.008):
    """Four hosts with wells and small discs in a thermal box."""
    ring = RingGeometry(.12, .15, .05, well_radius=.09, well_depth=depth)
    rng = np.random.default_rng(seed)
    centres = [(.4, .3), (1., .7), (1.6, .3), (1.55, .75)]
    specs = [host(*p, *rng.normal(0, math.sqrt(temperature/25), 2),
                  angle=rng.uniform(-3, 3), ring=ring) for p in centres]
    while len(specs) < 4 + discs:
        p = rng.uniform(radius + .01, [2 - radius - .01, 1 - radius - .01])
        if any(np.hypot(*(p - np.array(q))) < .15 + radius + .01 for q in centres):
            continue
        if any(np.hypot(*(p - np.array(s.position))) < 2*radius + .002 for s in specs[4:]):
            continue
        specs.append(disc(*p, *rng.normal(0, math.sqrt(temperature), 2), r=radius))
    sim = make(specs, {i: ring for i in range(4)}, temperature=temperature)
    return sim, ring


def mean_occupancy(sim, start, samples, step=.5):
    sim.advance_to(start)
    total = 0.
    for k in range(1, samples + 1):
        sim.advance_to(start + step*k)
        total += sim.occupancy().mean()
    return total/samples


def ideal_mean_occupancy(ring, discs, temperature, radius, width=2., height=1., hosts=4):
    """Ideal points: each disc independently in a given well with probability p."""
    well = math.pi*ring.well_radius**2
    cavity = math.pi*(ring.inner_radius - radius)**2
    free = ((width - 2*radius)*(height - 2*radius) - hosts*math.pi*(ring.outer_radius + radius)**2
            + hosts*(cavity - well))
    weight = well*math.exp(ring.well_depth/temperature)
    return discs*weight/(free + hosts*weight)


@unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
class StepLawTests(unittest.TestCase):
    def test_entering_the_well_follows_the_step_law(self):
        # Mouth faces the disc; it runs radially through the mouth to the well.
        sim = make([host(1., .5, angle=math.pi), disc(.5, .5, 1.)], {0: WELL})
        e0, p0 = sim.energy(), sim.world.bodies.momentum().copy()
        l0 = sim.world.bodies.angular_momentum()
        sim.step_collision()
        event = sim.events[-1]
        self.assertEqual((event.kind, event.metadata["step"]), ("step", "in"))
        self.assertAlmostEqual(sim.time, (1. - .08) - .5, delta=1e-14)
        self.assertEqual(list(sim.occupancy()), [1])
        # g' = -sqrt(g^2 + 2 D eps), D = 1 + 1/25, with g = -1 (inward).
        d = 1. + 1/25.
        g_out = -math.sqrt(1. + 2*d*1.5)
        v = sim.world.bodies.vel
        self.assertAlmostEqual(v[1, 0] - v[0, 0], -g_out, delta=1e-14)
        self.assertAlmostEqual(sim.energy(), e0, delta=1e-14)
        np.testing.assert_allclose(sim.world.bodies.momentum(), p0, atol=1e-14)
        self.assertAlmostEqual(sim.world.bodies.angular_momentum(), l0, delta=1e-14)
        self.assertEqual(sim.world.bodies.omega[0], 0.)     # radial: no torque
        np.testing.assert_array_equal(sim.host_statistics(), [[1, 0, 0]])

    def test_leaving_needs_enough_radial_energy(self):
        # From the centre outward at speed u: leaves iff u^2 > 2 D eps (D = 1 + 1/25).
        threshold = math.sqrt(2*(1 + 1/25.)*1.5)
        for u, crossed in ((.9*threshold, False), (1.1*threshold, True)):
            with self.subTest(speed=u):
                sim = make([host(1., .5, angle=0.), disc(1., .5, u)], {0: WELL})
                self.assertEqual(list(sim.occupancy()), [1])
                e0 = sim.energy()
                sim.step_collision()
                event = sim.events[-1]
                self.assertEqual(event.kind, "step")
                self.assertEqual(event.metadata["step"], "out" if crossed else "refused")
                self.assertEqual(list(sim.occupancy()), [0 if crossed else 1])
                self.assertAlmostEqual(sim.energy(), e0, delta=1e-13)
                rel = sim.world.bodies.vel[1, 0] - sim.world.bodies.vel[0, 0]
                expected = (math.sqrt(u*u - 2*(1 + 1/25.)*1.5) if crossed else -u)
                self.assertAlmostEqual(rel, expected, delta=1e-13)

    def test_zero_depth_well_only_counts(self):
        ring = RingGeometry(.12, .15, .08, well_radius=.08)
        sim = make([host(1., .5, angle=math.pi), disc(.5, .5, 1.)], {0: ring})
        sim.step_collision()
        self.assertEqual(sim.events[-1].metadata["step"], "in")
        np.testing.assert_allclose(sim.world.bodies.vel[1], [1., 0.])
        np.testing.assert_allclose(sim.world.bodies.vel[0], [0., 0.])

    def test_well_must_clear_the_wall(self):
        ring = RingGeometry(.12, .15, .08, well_radius=.11)
        world = World(BodyArrays.from_specs([host(1., .5, ring=ring), disc(.3, .5)]),
                      box_walls(2., 1.), rings={0: ring})
        self.assertIn("clear of the wall", edmd_unsupported_reason(world))
        with self.assertRaises(ValueError):
            RingGeometry(.12, .15, .08, well_depth=1.)


@unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
class WellDynamicsTests(unittest.TestCase):
    def test_first_law_with_wells_over_a_long_thermal_run(self):
        sim, _ = occupancy_box(1.5, 1.)
        for k in range(1, 201):
            snap = sim.advance_to(k*.5)
        stats = sim.host_statistics()
        self.assertGreater(stats[:, 0].sum(), 1000)
        self.assertGreater(stats[:, 2].sum(), 100)        # refused exits happen
        self.assertLessEqual(abs(snap.energy_residual), 1e-10*abs(snap.energy))
        # All discs start outside: entries minus exits, from the log and from
        # the counters, equal the discs bound now.
        net = sum(1 if e.metadata["step"] == "in" else -1 for e in sim.events
                  if e.kind == "step" and e.metadata["step"] != "refused")
        self.assertEqual(net, sim.occupancy().sum())
        self.assertEqual(stats[:, 0].sum() - stats[:, 1].sum(), sim.occupancy().sum())

    def test_a_frozen_smooth_host_traps_forever_but_a_mobile_one_releases(self):
        # One disc in the well whose radial energy at the step is too small.
        def trapped(mass, inertia, roughness, temperature):
            specs = [host(1., .5, angle=0., m=mass, inertia=inertia), disc(1.03, .52, .6, 1.1)]
            return make(specs, {0: WELL}, roughness=roughness, temperature=temperature)
        frozen = trapped(1e12, 1e12, 0., None)
        frozen.advance_to(2000.)
        stats = frozen.host_statistics()
        self.assertEqual(stats[0, 1], 0)                  # never left
        self.assertGreater(stats[0, 2], 100)              # but tried many times
        mobile = trapped(25., None, 1., 1.)
        mobile.advance_to(2000.)
        self.assertGreater(mobile.host_statistics()[0, 1], 0)

    def test_hard_only_occupancy_matches_reference_and_ignores_temperature(self):
        means = {}
        for temperature in (.5, 2.):
            sim, ring = occupancy_box(0., temperature)
            means[temperature] = mean_occupancy(sim, 50., 1600)
        reference = ideal_mean_occupancy(ring, 60, 1., .008)
        for temperature, mean in means.items():
            with self.subTest(temperature=temperature):
                self.assertAlmostEqual(mean, reference, delta=.1*reference)
        self.assertAlmostEqual(means[.5], means[2.], delta=.12*reference)

    def test_binding_raises_occupancy_as_the_reference_predicts(self):
        sim, ring = occupancy_box(1.5, 1.)
        mean = mean_occupancy(sim, 50., 1600)
        reference = ideal_mean_occupancy(ring, 60, 1., .008)
        self.assertGreater(reference, 3*ideal_mean_occupancy(
            RingGeometry(.12, .15, .05, well_radius=.09), 60, 1., .008))
        # Finite disc size crowds the small well: a few per cent below ideal points.
        self.assertAlmostEqual(mean, reference, delta=.12*reference)
        self.assertLess(mean, reference)


if __name__ == "__main__":
    unittest.main()
