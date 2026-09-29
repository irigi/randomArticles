import math
import unittest

import numpy as np

from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.core.broadphase import swept_pairs
from microthermo.core.numeric import numba_available, polygon_witnesses, reach_mask
from microthermo.core.geometry import convex_witnesses, world_polygon
from microthermo.core.ccd import body_pair_toi
from microthermo.core.state import BodyArrays, BodySpec, Shape
from microthermo.core.boundaries import Wall
from microthermo.runner.simulation import Simulation, World


class BroadphaseTests(unittest.TestCase):
    def test_swept_grid_contains_every_circumcircle_candidate(self):
        rng = np.random.default_rng(441)
        for count in (32, 128):
            pos = rng.uniform(-2, 2, (count, 2))
            vel = rng.normal(0, 3, (count, 2))
            radius = rng.uniform(.001, .08, count)
            horizon = .07
            candidates = set(swept_pairs(pos, vel, radius, horizon))
            for a in range(count):
                for b in range(a + 1, count):
                    dp, dv = pos[b] - pos[a], vel[b] - vel[a]
                    speed_sq = float(dv @ dv)
                    closest_time = np.clip(-float(dp @ dv) / speed_sq, 0, horizon) if speed_sq else 0
                    if np.linalg.norm(dp + closest_time * dv) <= radius[a] + radius[b]:
                        self.assertIn((a, b), candidates)

    def test_area_fraction_is_constant_and_dilute(self):
        for preset, area_factor in (("carnot_triangles", 3 * math.sqrt(3) / 4),
                                    ("carnot_discs", math.pi)):
            fractions = []
            for count in (16, 48, 128):
                sim = load_preset(RunConfig(preset=preset, particles=count))
                radius = sim.world.bodies.radius[0]
                fractions.append(count * area_factor * radius**2 / 1.25)
                self.assertAlmostEqual(sim.world.metadata["particle_area_fraction_min"],
                                       fractions[-1])
            np.testing.assert_allclose(fractions, fractions[0], rtol=1e-14)
            self.assertLess(fractions[0], .1)

    def test_grid_and_all_pairs_produce_same_carnot_trajectory(self):
        for preset in ("carnot_discs", "carnot_triangles"):
            config = RunConfig(preset=preset, particles=32, seed=123, max_horizon=.02)
            grid = load_preset(config)
            reference = load_preset(config)
            grid.numeric_backend = "python"
            reference.pair_search = "all"
            reference.numeric_backend = "python"
            a = grid.advance_to(.15)
            b = reference.advance_to(.15)
            self.assertEqual([(e.time, e.kind, e.participants) for e in grid.events],
                             [(e.time, e.kind, e.participants) for e in reference.events])
            np.testing.assert_array_equal(a.position, b.position)
            np.testing.assert_array_equal(a.velocity, b.velocity)
            self.assertEqual(grid.ledger.heat_hot.value, reference.ledger.heat_hot.value)
            self.assertEqual(grid.ledger.work_on.value, reference.ledger.work_on.value)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_numba_reach_mask_matches_python(self):
        rng = np.random.default_rng(331)
        pos = rng.normal(size=(64, 2))
        vel = rng.normal(size=(64, 2))
        radius = rng.uniform(.01, .05, 64)
        omega = rng.normal(size=64)
        pairs = np.asarray(swept_pairs(pos, vel, radius, .03), dtype=np.int64)
        np.testing.assert_array_equal(
            reach_mask(pos, vel, radius, omega, pairs, .03, "python"),
            reach_mask(pos, vel, radius, omega, pairs, .03, "numba"))

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_numba_triangle_witnesses_match_reference(self):
        rng = np.random.default_rng(71)
        local = np.array([[0., 1.], [-math.sqrt(3)/2, -.5],
                          [math.sqrt(3)/2, -.5]]) * .07
        for _ in range(200):
            a = world_polygon(local, rng.normal(0, .12, 2), rng.uniform(-math.pi, math.pi))
            b = world_polygon(local, rng.normal(0, .12, 2), rng.uniform(-math.pi, math.pi))
            reference = convex_witnesses(a, b)
            compiled = polygon_witnesses(a, b, "numba")
            self.assertAlmostEqual(reference[0], compiled[0], delta=2e-15)
            np.testing.assert_allclose(reference[1], compiled[1], atol=2e-14, rtol=0)
            np.testing.assert_allclose(reference[2], compiled[2], atol=2e-14, rtol=0)
            np.testing.assert_allclose(reference[3], compiled[3], atol=2e-14, rtol=0)
            self.assertEqual(reference[4:], compiled[4:])

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_compiled_triangle_trajectory_matches_reference(self):
        config = RunConfig(preset="carnot_triangles", particles=32, seed=123,
                           max_horizon=.02, numeric_backend="numba")
        compiled = load_preset(config)
        reference = load_preset(config)
        reference.numeric_backend = "python"
        a = compiled.advance_to(.15)
        b = reference.advance_to(.15)
        self.assertEqual([(e.kind, e.participants) for e in compiled.events],
                         [(e.kind, e.participants) for e in reference.events])
        np.testing.assert_allclose([e.time for e in compiled.events],
                                   [e.time for e in reference.events], atol=2e-12, rtol=0)
        np.testing.assert_allclose(a.position, b.position, atol=2e-11, rtol=0)
        np.testing.assert_allclose(a.velocity, b.velocity, atol=2e-11, rtol=0)
        self.assertAlmostEqual(compiled.ledger.heat_hot.value,
                               reference.ledger.heat_hot.value, delta=2e-11)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_compiled_polygon_ccd_matches_reference_cases(self):
        rng = np.random.default_rng(917)
        cases = [
            (0.1, 0.0, 0.0, 0.0, 0.0, 0.0),  # touching and stationary
            (0.1, 0.0, 0.0, 0.0, 9.0, -9.0),  # touching and rotating
            (0.25, 0.0, 1000.0, -1000.0, 0.0, 0.0),
        ]
        for _ in range(30):
            cases.append((rng.uniform(.12, .45), rng.uniform(-.06, .06),
                          rng.uniform(-4, 4), rng.uniform(-4, 4),
                          rng.uniform(-20, 20), rng.uniform(-20, 20)))
        for distance, lateral, va, vb, wa, wb in cases:
            state = BodyArrays.from_specs([
                BodySpec((0., 0.), (va, 0.), radius=.05, shape=Shape.TRIANGLE,
                         angle=.2, omega=wa),
                BodySpec((distance, lateral), (vb, 0.), radius=.05,
                         shape=Shape.TRIANGLE, angle=-.3, omega=wb),
            ])
            python = body_pair_toi(state, 0, 1, .05, numeric_backend="python")
            compiled = body_pair_toi(state, 0, 1, .05, numeric_backend="numba")
            self.assertEqual(python.status, compiled.status)
            self.assertEqual(python.reason, compiled.reason)
            if python.contact is not None:
                self.assertAlmostEqual(python.time, compiled.time, delta=2e-11)
                np.testing.assert_allclose(python.contact.point, compiled.contact.point,
                                           atol=2e-10, rtol=0)
                np.testing.assert_allclose(python.contact.normal, compiled.contact.normal,
                                           atol=2e-10, rtol=0)
                self.assertEqual((python.contact.feature_a, python.contact.feature_b),
                                 (compiled.contact.feature_a, compiled.contact.feature_b))

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_compiled_checkpoint_and_sampling_keep_event_chronology(self):
        config = RunConfig(preset="carnot_triangles", particles=32, seed=124,
                           max_horizon=.02, numeric_backend="numba")
        direct = load_preset(config)
        direct_end = direct.advance_to(.15)
        sampled = load_preset(config)
        sampled.advance_to(.037)
        checkpoint = sampled.checkpoint()
        for t in (.051, .079, .113, .15):
            sampled_end = sampled.advance_to(t)
        events = [(e.time, e.kind, e.participants) for e in sampled.events]
        self.assertEqual([(e.time, e.kind, e.participants) for e in direct.events],
                         events)
        np.testing.assert_array_equal(direct_end.position, sampled_end.position)
        sampled.restore(checkpoint)
        replay_end = sampled.advance_to(.15)
        self.assertEqual(events,
                         [(e.time, e.kind, e.participants) for e in sampled.events])
        np.testing.assert_array_equal(sampled_end.position, replay_end.position)
        np.testing.assert_array_equal(sampled_end.velocity, replay_end.velocity)
        self.assertLessEqual(sampled.max_penetration, sampled.tol.geometry)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_bounded_walls_match_unpruned_reference(self):
        config = RunConfig(preset="carnot_triangles", particles=32, seed=124,
                           max_horizon=.02, numeric_backend="numba")
        bounded = load_preset(config)
        reference = load_preset(config)
        reference.wall_search = "all"
        a = bounded.advance_to(.15)
        b = reference.advance_to(.15)
        self.assertEqual([(e.time, e.kind, e.participants) for e in bounded.events],
                         [(e.time, e.kind, e.participants) for e in reference.events])
        np.testing.assert_array_equal(a.position, b.position)
        np.testing.assert_array_equal(a.velocity, b.velocity)
        self.assertEqual(bounded.ledger.heat_hot.value, reference.ledger.heat_hot.value)
        self.assertEqual(bounded.max_penetration, reference.max_penetration)

    def test_bounded_wall_checks_keep_near_and_fast_contacts(self):
        config = RunConfig(preset="carnot_triangles", particles=1, seed=123)
        bounded = load_preset(config)
        reference = load_preset(config)
        wall = bounded.world.walls[3]  # moving piston
        for sim in (bounded, reference):
            state = sim.world.bodies
            state.pos[0] = (wall.point_at(0.)[0] - state.radius[0] - .01, .5)
            state.vel[0] = (1000., 0.)
        reference.wall_search = "all"
        speed = wall.speed_bound(0., .001)
        start = wall.point_at(0.)
        fast = bounded._wall_toi(0, 3, .001, speed, start)
        direct = reference._wall_toi(0, 3, .001)
        self.assertEqual(fast.status, direct.status)
        self.assertAlmostEqual(fast.time, direct.time, delta=1e-12)

        # A penetrated wall must remain visible to the overlap diagnostic.
        for sim in (bounded, reference):
            sim.world.bodies.pos[0, 0] = .01
        self.assertEqual(bounded._penetration()[1], reference._penetration()[1])
        self.assertAlmostEqual(bounded._penetration()[0],
                               reference._penetration()[0], delta=1e-12)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_fixed_wall_kernel_matches_reference_contacts(self):
        for x, vx, omega, wall_speed, horizon in (
            (.3, -1000., 0., 0., .001),
            (.08, 0., 20., 0., .1),
            (.3, -1., 3., .2, .5),
        ):
            with self.subTest(vx=vx, omega=omega, wall_speed=wall_speed):
                specs = [BodySpec((x, .5), (vx, 0.), radius=.05,
                                  shape=Shape.TRIANGLE, angle=.2, omega=omega)]
                wall = Wall(np.array([0., .5]), np.array([1., 0.]), 1.,
                            velocity=np.array([wall_speed, 0.]))
                compiled = Simulation(World(BodyArrays.from_specs(specs), [wall]),
                                      numeric_backend="numba", wall_search="all")
                reference = Simulation(World(BodyArrays.from_specs(specs), [wall]),
                                       numeric_backend="python", wall_search="all")
                a = compiled._wall_toi(0, 0, horizon)
                b = reference._wall_toi(0, 0, horizon)
                self.assertEqual(a.status, b.status)
                self.assertEqual(a.reason, b.reason)
                if a.contact is not None:
                    self.assertAlmostEqual(a.time, b.time, delta=2e-11)
                    np.testing.assert_allclose(a.contact.point, b.contact.point,
                                               atol=2e-10, rtol=0)
                    self.assertEqual(a.contact.feature_a, b.contact.feature_a)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_fixed_wall_kernel_matches_carnot_trajectory(self):
        config = RunConfig(preset="carnot_triangles", particles=48, seed=123,
                           max_horizon=.02, numeric_backend="numba",
                           wall_search="bounded")
        compiled = load_preset(config)
        reference = load_preset(config)
        reference.wall_kernel = "python"
        a = compiled.advance_to(.12)
        b = reference.advance_to(.12)
        self.assertEqual([(e.kind, e.participants) for e in compiled.events],
                         [(e.kind, e.participants) for e in reference.events])
        np.testing.assert_allclose([e.time for e in compiled.events],
                                   [e.time for e in reference.events], atol=2e-12, rtol=0)
        np.testing.assert_allclose(a.position, b.position, atol=2e-11, rtol=0)
        np.testing.assert_allclose(a.velocity, b.velocity, atol=2e-11, rtol=0)
        self.assertAlmostEqual(compiled.ledger.heat_hot.value,
                               reference.ledger.heat_hot.value, delta=2e-11)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_compiled_penetration_bounds_match_python_near_all_carnot_walls(self):
        config = RunConfig(preset="carnot_triangles", particles=1,
                           numeric_backend="numba", wall_search="bounded")
        compiled = load_preset(config)
        reference = load_preset(config)
        reference.penetration_kernel = "python"
        piston_x = compiled.world.walls[3].point_at(0.)[0]
        for position in ((.01, .5), (.5, .01), (.5, .99),
                         (piston_x - .01, .5), (.5, .5)):
            for sim in (compiled, reference):
                sim.world.bodies.pos[0] = position
            a = compiled._penetration()
            b = reference._penetration()
            self.assertEqual(a[1], b[1])
            self.assertAlmostEqual(a[0], b[0], delta=2e-14)


if __name__ == "__main__":
    unittest.main()
