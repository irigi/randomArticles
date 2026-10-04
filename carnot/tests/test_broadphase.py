import math
import unittest

import numpy as np

from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.core.broadphase import (swept_pairs, swept_pairs_sweep,
                                         swept_pairs_sweep_array)
from microthermo.core.numeric import (controlled_cam_disc_toi,
                                      controlled_cam_triangle_toi, disc_pairs_toi,
                                      fixed_disc_walls_toi,
                                      _polygon_gap_at_numba,
                                      _polygon_witness_at_numba,
                                      numba_available, polygon_pair_toi,
                                      polygon_pairs_toi, polygon_witnesses,
                                      reach_mask, swept_circle_mask,
                                      triangle_pair_gaps,
                                      wall_reach_mask)
from microthermo.core.geometry import convex_separation, convex_witnesses, world_polygon
from microthermo.core.ccd import body_pair_toi
from microthermo.core.events import TOIStatus
from microthermo.core.state import BodyArrays, BodySpec, Shape, Tolerances
from microthermo.core.boundaries import Wall
from microthermo.runner.simulation import Simulation, World
from microthermo.measurements.ledger import EnergyLedger


class BroadphaseTests(unittest.TestCase):
    def test_packed_triangle_geometry_survives_independent_checkpoint_copy(self):
        state = load_preset(RunConfig(preset="carnot_triangles", particles=32,
                                      seed=123)).world.bodies
        self.assertIsNotNone(state.packed_triangles)
        saved = state.copy()
        original_vertex = float(saved.polygons[0][0, 0])
        state.polygons[0][0, 0] += .01
        self.assertAlmostEqual(state.packed_triangles[0, 0, 0],
                               original_vertex + .01)
        self.assertAlmostEqual(saved.polygons[0][0, 0], original_vertex)
        self.assertAlmostEqual(saved.packed_triangles[0, 0, 0], original_vertex)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_polygon_sat_only_gap_matches_full_witness(self):
        state = load_preset(RunConfig(preset="carnot_triangles", particles=2,
                                      numeric_backend="numba")).world.bodies
        rng = np.random.default_rng(833)
        for _ in range(100):
            state.pos[0] = rng.normal(size=2)
            state.pos[1] = state.pos[0] + rng.normal(0., .08, size=2)
            state.vel[:] = rng.normal(size=(2, 2))
            state.angle[:] = rng.uniform(-math.pi, math.pi, size=2)
            state.omega[:] = rng.normal(0., 5., size=2)
            t = rng.uniform(0., .1)
            args = (state.polygons[0], state.polygons[1],
                    state.pos[0], state.pos[1], state.vel[0], state.vel[1],
                    state.angle[0], state.angle[1], state.omega[0],
                    state.omega[1], t)
            self.assertAlmostEqual(_polygon_gap_at_numba(*args),
                                   _polygon_witness_at_numba(*args)[0],
                                   delta=2e-14)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_polygon_batch_decoding_keeps_simultaneous_and_failure_order(self):
        specs = [BodySpec((0., 0.), (1., 0.), radius=.1, shape=Shape.TRIANGLE),
                 BodySpec((1., 0.), (0., 0.), radius=.1, shape=Shape.TRIANGLE),
                 BodySpec((0., 2.), (1., 0.), radius=.1, shape=Shape.TRIANGLE),
                 BodySpec((1., 2.), (0., 0.), radius=.1, shape=Shape.TRIANGLE)]
        fast = Simulation(World(BodyArrays.from_specs(specs), []),
                          numeric_backend="numba", pair_search="all")
        reference = Simulation(World(BodyArrays.from_specs(specs), []),
                               numeric_backend="numba", pair_search="all",
                               pair_kernel="scalar")
        a, b = fast._earliest(1.), reference._earliest(1.)
        self.assertEqual(a.status, b.status)
        self.assertAlmostEqual(a.time, b.time, delta=2e-12)
        self.assertEqual([(c.a, c.b) for c in a.contacts],
                         [(c.a, c.b) for c in b.contacts])
        self.assertEqual([(c.a, c.b) for c in a.contacts], [(0, 1), (2, 3)])
        for sim in (fast, reference):
            sim.world.bodies.pos[1] = sim.world.bodies.pos[0]
        failure, expected = fast._earliest(1.), reference._earliest(1.)
        self.assertEqual(failure.status, expected.status)
        self.assertEqual(failure.reason, expected.reason)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_swept_circle_filter_is_conservative_and_matches_python(self):
        rng = np.random.default_rng(1729)
        for horizon in (0., .03, .2):
            pos = rng.normal(size=(80, 2))
            vel = rng.normal(0., 12., size=(80, 2))
            radius = rng.uniform(.002, .08, size=80)
            pairs = np.array([(a, b) for a in range(80)
                              for b in range(a + 1, 80)], dtype=np.int64)
            fast = swept_circle_mask(pos, vel, radius, pairs, horizon,
                                     1e-10, "numba")
            slow = swept_circle_mask(pos, vel, radius, pairs, horizon,
                                     1e-10, "python")
            np.testing.assert_array_equal(fast, slow)
            for t in (0., .25 * horizon, .5 * horizon, horizon):
                centers = pos + t * vel
                delta = centers[pairs[:, 1]] - centers[pairs[:, 0]]
                squared = np.sum(delta * delta, axis=1)
                limit = (radius[pairs[:, 0]] + radius[pairs[:, 1]])**2
                self.assertTrue(np.all(fast[squared <= limit]))
        pos = np.array([[0., 0.], [1., .1], [2., 0.]])
        vel = np.array([[0., 0.], [-2., 0.], [-1000., 0.]])
        radius = np.array([.05, .05, .05])
        pairs = np.array([[0, 1], [0, 2]], dtype=np.int64)
        self.assertTrue(np.all(swept_circle_mask(pos, vel, radius, pairs,
                                                  1., 1e-10, "numba")))

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_swept_circle_filter_preserves_triangle_events(self):
        config = RunConfig(preset="carnot_triangles", particles=96, seed=128,
                           pair_search="sweep", numeric_backend="numba")
        filtered = load_preset(config)
        reference = load_preset(config)
        reference.pair_kernel = "scalar"
        a = filtered.advance_to(.2)
        b = reference.advance_to(.2)
        self.assertEqual([(e.kind, e.participants) for e in filtered.events],
                         [(e.kind, e.participants) for e in reference.events])
        np.testing.assert_allclose([e.time for e in filtered.events],
                                   [e.time for e in reference.events], atol=2e-11, rtol=0)
        np.testing.assert_allclose(a.position, b.position, atol=2e-9, rtol=0)
        np.testing.assert_allclose(a.velocity, b.velocity, atol=2e-9, rtol=0)
        self.assertAlmostEqual(filtered.ledger.heat_hot.value,
                               reference.ledger.heat_hot.value, delta=2e-9)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_controlled_cam_triangle_batch_matches_scalar_queries(self):
        sim = load_preset(RunConfig(preset="carnot_triangles", particles=24,
                                    numeric_backend="numba"))
        rng = np.random.default_rng(9217)
        wall = sim.world.walls[3]
        state = sim.world.bodies
        for phase in (.2, 1.4, 2.0, 4.6):
            sim.world.mechanism.shaft.phi = phase
            piston_x = wall.point_at(sim.time)[0]
            for a in range(state.n):
                state.pos[a] = (piston_x - state.radius[a] - rng.uniform(.001, .2),
                                rng.uniform(.1, .9))
                state.vel[a] = (rng.uniform(-4, 8), rng.uniform(-3, 3))
                state.angle[a] = rng.uniform(-math.pi, math.pi)
                state.omega[a] = rng.uniform(-10, 10)
            horizon = .03
            speed = wall.speed_bound(sim.time, sim.time + horizon)
            batch = controlled_cam_triangle_toi(
                state, np.ones(state.n, dtype=np.bool_),
                sim.world.mechanism, speed, horizon, sim.tol)
            for a in range(state.n):
                scalar = sim._wall_toi(a, 3, horizon, speed, wall.point_at(sim.time))
                self.assertEqual(int(batch[a, 0]), {"no_collision": 0,
                                                    "collision": 1,
                                                    "indeterminate": 2}[scalar.status.value])
                if scalar.contact:
                    self.assertAlmostEqual(batch[a, 1], scalar.time, delta=2e-11)
                    np.testing.assert_allclose(batch[a, 3:5], scalar.contact.point,
                                               atol=2e-10, rtol=0)
                    self.assertEqual(int(batch[a, 5])*2, scalar.contact.feature_a)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_triangle_cam_contact_across_each_branch_boundary(self):
        sim = load_preset(RunConfig(preset="carnot_triangles", particles=16,
                                    shaft_speed=1.5, numeric_backend="numba"))
        state, wall = sim.world.bodies, sim.world.walls[3]
        horizon = .2
        for boundary, gap in ((math.pi/2,.08), (math.pi,.08),
                              (3*math.pi/2,.15), (2*math.pi,.15)):
            with self.subTest(boundary=boundary):
                sim.world.mechanism.shaft.phi = boundary-.12
                state.pos[0] = (wall.point_at(sim.time)[0]-state.radius[0]-gap,.5)
                state.vel[0] = (1.,0.)
                state.angle[0] = .2
                state.omega[0] = 1.
                wall_start = wall.point_at(sim.time)
                speed = wall.speed_bound(sim.time,sim.time+horizon)
                compiled = controlled_cam_triangle_toi(
                    state, np.arange(state.n)==0, sim.world.mechanism,
                    speed, horizon, sim.tol)[0]
                bounded = sim._wall_toi(0,3,horizon,speed,wall_start)
                unpruned = sim._wall_toi(0,3,horizon)
                self.assertEqual(int(compiled[0]),1)
                self.assertEqual(bounded.status,TOIStatus.COLLISION)
                self.assertEqual(unpruned.status,TOIStatus.COLLISION)
                self.assertGreater(compiled[1],.12/1.5)
                self.assertAlmostEqual(compiled[1],bounded.time,delta=2e-11)
                self.assertAlmostEqual(compiled[1],unpruned.time,delta=2e-11)
                np.testing.assert_allclose(compiled[3:5],bounded.contact.point,
                                           atol=2e-10,rtol=0)
                self.assertEqual(int(compiled[5])*2,bounded.contact.feature_a)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_controlled_cam_triangle_batch_preserves_trajectory(self):
        config = RunConfig(preset="carnot_triangles", particles=96, seed=125,
                           pair_search="sweep", numeric_backend="numba")
        compiled = load_preset(config)
        reference = load_preset(config)
        reference.cam_kernel = "python"
        a = compiled.advance_to(.5)
        b = reference.advance_to(.5)
        self.assertGreater(sum(bool(e.metadata and e.metadata.get("boundary") == 3)
                               for e in compiled.events), 0)
        self.assertEqual([(e.kind, e.participants) for e in compiled.events],
                         [(e.kind, e.participants) for e in reference.events])
        np.testing.assert_allclose([e.time for e in compiled.events],
                                   [e.time for e in reference.events], atol=2e-11, rtol=0)
        np.testing.assert_allclose(a.position, b.position, atol=2e-9, rtol=0)
        np.testing.assert_allclose(a.velocity, b.velocity, atol=2e-9, rtol=0)
        self.assertAlmostEqual(compiled.ledger.heat_hot.value,
                               reference.ledger.heat_hot.value, delta=2e-9)
        self.assertAlmostEqual(compiled.ledger.work_on.value,
                               reference.ledger.work_on.value, delta=2e-9)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_triangle_pair_gap_batch_matches_reference(self):
        sim = load_preset(RunConfig(preset="carnot_triangles", particles=32,
                                    numeric_backend="numba", pair_search="sweep"))
        s = sim.world.bodies
        for offset in (0., .5 * s.radius[0], 2.5 * s.radius[0]):
            s.pos[1] = s.pos[0] + (offset, 0.)
            pairs = np.array([[0, 1], [0, 2]], dtype=np.int64)
            gaps = triangle_pair_gaps(s, pairs)
            for i, (a, b) in enumerate(pairs):
                if np.linalg.norm(s.pos[b] - s.pos[a]) > s.radius[a] + s.radius[b]:
                    self.assertTrue(np.isinf(gaps[i]))
                else:
                    pa = world_polygon(s.polygons[a], s.pos[a], s.angle[a])
                    pb = world_polygon(s.polygons[b], s.pos[b], s.angle[b])
                    self.assertAlmostEqual(gaps[i], convex_separation(pa, pb)[0],
                                           delta=2e-12)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_triangle_penetration_batch_reports_same_overlap(self):
        config = RunConfig(preset="carnot_triangles", particles=32, seed=29,
                           pair_search="sweep", numeric_backend="numba")
        compiled = load_preset(config)
        reference = load_preset(config)
        reference.penetration_kernel = "python"
        for sim in (compiled, reference):
            sim.world.bodies.pos[1] = sim.world.bodies.pos[0]
        a, b = compiled._penetration(), reference._penetration()
        self.assertEqual(a[1], b[1])
        self.assertLess(a[0], 0.)
        self.assertAlmostEqual(a[0], b[0], delta=2e-12)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_triangle_penetration_batch_preserves_longer_trajectory(self):
        config = RunConfig(preset="carnot_triangles", particles=96, seed=127,
                           pair_search="sweep", numeric_backend="numba")
        compiled = load_preset(config)
        reference = load_preset(config)
        reference.penetration_kernel = "python"
        a = compiled.advance_to(.2)
        b = reference.advance_to(.2)
        self.assertEqual([(e.kind, e.participants) for e in compiled.events],
                         [(e.kind, e.participants) for e in reference.events])
        np.testing.assert_allclose([e.time for e in compiled.events],
                                   [e.time for e in reference.events], atol=2e-12, rtol=0)
        np.testing.assert_allclose(a.position, b.position, atol=2e-10, rtol=0)
        np.testing.assert_allclose(a.velocity, b.velocity, atol=2e-10, rtol=0)
        self.assertAlmostEqual(a.energy_residual, b.energy_residual, delta=2e-10)
        self.assertAlmostEqual(compiled.max_penetration,
                               reference.max_penetration, delta=2e-12)

    def test_checkpoint_ledger_copy_preserves_counter_corrections(self):
        ledger = EnergyLedger(3.)
        ledger.heat_hot.value = 1.25
        ledger.heat_hot.correction = -2e-16
        ledger.support_impulse_y.value = .5
        clone = ledger.copy()
        self.assertEqual(clone.heat_hot.value, ledger.heat_hot.value)
        self.assertEqual(clone.heat_hot.correction, ledger.heat_hot.correction)
        self.assertEqual(clone.support_impulse_y.value, ledger.support_impulse_y.value)
        clone.heat_hot.add(.125)
        self.assertEqual(ledger.heat_hot.value, 1.25)

    def test_carnot_checkpoint_shares_frozen_cam_but_restores_shaft(self):
        sim = load_preset(RunConfig(preset="carnot_discs", particles=16))
        original_phase = sim.world.mechanism.shaft.phi
        self.assertFalse(sim.world.mechanism.cam.boundaries.flags.writeable)
        self.assertEqual(sim.energy(), sim.world.energy())
        checkpoint = sim.checkpoint()
        self.assertIs(checkpoint.mechanism.cam, sim.world.mechanism.cam)
        self.assertIsNot(checkpoint.mechanism.shaft, sim.world.mechanism.shaft)
        sim.world.mechanism.shaft.phi += .2
        sim.restore(checkpoint)
        self.assertEqual(sim.world.mechanism.shaft.phi, original_phase)
        self.assertIs(sim.world.mechanism.cam, checkpoint.mechanism.cam)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_disc_batch_decoding_keeps_simultaneous_and_failure_order(self):
        specs = [BodySpec((0., 0.), (1., 0.), radius=.1),
                 BodySpec((1., 0.), (0., 0.), radius=.1),
                 BodySpec((0., 2.), (1., 0.), radius=.1),
                 BodySpec((1., 2.), (0., 0.), radius=.1)]
        fast = Simulation(World(BodyArrays.from_specs(specs), []),
                          numeric_backend="numba", pair_search="all")
        reference = Simulation(World(BodyArrays.from_specs(specs), []),
                               numeric_backend="numba", pair_search="all",
                               pair_kernel="scalar")
        a, b = fast._earliest(1.), reference._earliest(1.)
        self.assertEqual(a.status, b.status)
        self.assertEqual(a.time, b.time)
        self.assertEqual(a.contacts, b.contacts)
        fast.world.bodies.pos[1] = fast.world.bodies.pos[0] + (1e-12, 0.)
        fast.world.bodies.vel[1] = (-100., 0.)
        failure = fast._earliest(1.)
        self.assertEqual(failure.status.value, "indeterminate")
        self.assertEqual(failure.reason, "pair (0, 1): coincident disc centers")

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_disc_batch_decoding_keeps_pair_before_wall_at_same_time(self):
        specs = [BodySpec((.1, .5), (-1., 0.), radius=.1),
                 BodySpec((.3, .5), (-2., 0.), radius=.1)]
        wall = Wall(np.array([0., .5]), np.array([1., 0.]), 1.)
        fast = Simulation(World(BodyArrays.from_specs(specs), [wall]),
                          numeric_backend="numba", pair_search="all")
        reference = Simulation(World(BodyArrays.from_specs(specs), [wall]),
                               numeric_backend="numba", pair_search="all",
                               pair_kernel="scalar")
        a, b = fast._earliest(.01), reference._earliest(.01)
        self.assertEqual(a.status, b.status)
        self.assertEqual(a.time, b.time)
        self.assertEqual(a.contacts, b.contacts)
        self.assertEqual([contact.boundary for contact in a.contacts], [None, 0])

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_disc_penetration_batch_matches_reference(self):
        config = RunConfig(preset="carnot_discs", particles=32, seed=391,
                           pair_search="sweep", numeric_backend="numba")
        compiled = load_preset(config)
        reference = load_preset(config)
        reference.penetration_kernel = "python"
        for mode in ("clear", "pair", "wall"):
            if mode == "pair":
                for sim in (compiled, reference):
                    sim.world.bodies.pos[1] = sim.world.bodies.pos[0] + (.01, 0.)
            elif mode == "wall":
                for sim in (compiled, reference):
                    sim.world.bodies.pos[1] = (0.001, .5)
            self.assertEqual(compiled._penetration()[1], reference._penetration()[1])
            self.assertAlmostEqual(compiled._penetration()[0],
                                   reference._penetration()[0], delta=2e-12)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_disc_penetration_batch_preserves_trajectory(self):
        config = RunConfig(preset="carnot_discs", particles=96, seed=126,
                           pair_search="sweep", numeric_backend="numba")
        compiled = load_preset(config)
        reference = load_preset(config)
        reference.penetration_kernel = "python"
        a = compiled.advance_to(.3)
        b = reference.advance_to(.3)
        self.assertEqual([(e.kind, e.participants) for e in compiled.events],
                         [(e.kind, e.participants) for e in reference.events])
        np.testing.assert_allclose([e.time for e in compiled.events],
                                   [e.time for e in reference.events], atol=2e-12, rtol=0)
        np.testing.assert_allclose(a.position, b.position, atol=2e-11, rtol=0)
        np.testing.assert_allclose(a.velocity, b.velocity, atol=2e-11, rtol=0)
        self.assertAlmostEqual(a.energy_residual, b.energy_residual, delta=2e-11)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_fixed_disc_wall_batch_matches_scalar_queries(self):
        cases = ((.3, -1000., 0., .001), (.08, 0., .2, .1),
                 (.3, -1., .2, .5), (.05, 1., 0., .1))
        for x, vx, wall_speed, horizon in cases:
            with self.subTest(x=x, vx=vx, wall_speed=wall_speed):
                state = BodyArrays.from_specs([BodySpec((x, .5), (vx, 0.), radius=.05)])
                wall = Wall(np.array([0., .5]), np.array([1., 0.]), 1.,
                            velocity=np.array([wall_speed, 0.]))
                sim = Simulation(World(state, [wall]), numeric_backend="numba",
                                 wall_search="all")
                raw = fixed_disc_walls_toi(
                    state, np.ones((1, 1), dtype=np.bool_),
                    np.array([wall.point_at(0.)]), np.array([wall.inward_normal]),
                    np.array([wall.velocity]), np.array([wall.speed_bound(0., horizon)]),
                    horizon, sim.tol)[0, 0]
                reference = sim._wall_toi(0, 0, horizon)
                self.assertEqual(int(raw[0]), {"no_collision": 0,
                                               "collision": 1,
                                               "indeterminate": 2}[reference.status.value])
                if reference.contact:
                    self.assertAlmostEqual(raw[1], reference.time, delta=2e-11)
                    np.testing.assert_allclose(raw[3:5], reference.contact.point,
                                               atol=2e-11, rtol=0)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_fixed_disc_wall_batch_preserves_carnot_trajectory(self):
        config = RunConfig(preset="carnot_discs", particles=96, seed=125,
                           pair_search="sweep", numeric_backend="numba")
        compiled = load_preset(config)
        reference = load_preset(config)
        reference.wall_kernel = "python"
        a = compiled.advance_to(.5)
        b = reference.advance_to(.5)
        self.assertEqual([(e.kind, e.participants) for e in compiled.events],
                         [(e.kind, e.participants) for e in reference.events])
        np.testing.assert_allclose([e.time for e in compiled.events],
                                   [e.time for e in reference.events], atol=2e-11, rtol=0)
        np.testing.assert_allclose(a.position, b.position, atol=2e-9, rtol=0)
        np.testing.assert_allclose(a.velocity, b.velocity, atol=2e-9, rtol=0)
        self.assertAlmostEqual(compiled.ledger.heat_hot.value,
                               reference.ledger.heat_hot.value, delta=2e-9)
        self.assertAlmostEqual(compiled.ledger.work_on.value,
                               reference.ledger.work_on.value, delta=2e-9)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_fixed_disc_wall_batch_excludes_finite_segments(self):
        config = RunConfig(preset="selective_membrane", particles=10,
                           duration=.2, max_horizon=.01, numeric_backend="numba")
        compiled = load_preset(config)
        reference = load_preset(config)
        reference.wall_kernel = "python"
        a = compiled.advance_to(.2)
        b = reference.advance_to(.2)
        self.assertEqual([(e.kind, e.participants) for e in compiled.events],
                         [(e.kind, e.participants) for e in reference.events])
        np.testing.assert_allclose(a.position, b.position, atol=2e-11, rtol=0)
        self.assertAlmostEqual(a.energy_residual, b.energy_residual, delta=2e-11)

    def test_triangle_default_count_doubles_without_overriding_explicit_counts(self):
        self.assertEqual(RunConfig(preset="carnot_triangles").particles, 96)
        self.assertEqual(RunConfig(preset="carnot_triangles", particles=32).particles, 32)
        self.assertEqual(RunConfig(preset="carnot_discs").particles, 32)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_batched_wall_reach_contains_scalar_bound(self):
        rng = np.random.default_rng(812)
        pos = rng.normal(size=(200, 2))
        vel = rng.normal(size=(200, 2))
        radius = rng.uniform(.002, .08, 200)
        points = rng.normal(size=(4, 2))
        normals = np.array([[1., 0.], [0., 1.], [-1., 0.], [0., -1.]])
        speeds = rng.uniform(0, 1, 4)
        segments = np.array([False, False, True, False])
        for horizon in (0., .01, .2):
            compiled = wall_reach_mask(pos, vel, radius, points, normals,
                                       speeds, segments, horizon, 1e-10, "numba")
            reference = wall_reach_mask(pos, vel, radius, points, normals,
                                        speeds, segments, horizon, 1e-10, "python")
            np.testing.assert_array_equal(compiled, reference)
            for a in range(len(radius)):
                for wi in range(len(speeds)):
                    if segments[wi]:
                        self.assertTrue(compiled[a, wi])
                        continue
                    gap = float((pos[a] - points[wi]) @ normals[wi] - radius[a])
                    speed = abs(float(vel[a] @ normals[wi]))
                    scalar = gap <= horizon * (speed + speeds[wi]) + 1e-10
                    if scalar:
                        self.assertTrue(compiled[a, wi])

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_batched_wall_reach_preserves_disc_trajectory(self):
        config = RunConfig(preset="carnot_discs", particles=96, seed=124,
                           pair_search="sweep", numeric_backend="numba")
        compiled = load_preset(config)
        reference = load_preset(config)
        compiled.pair_kernel = "scalar"
        reference.pair_kernel = "scalar"
        reference.numeric_backend = "python"
        a = compiled.advance_to(.05)
        b = reference.advance_to(.05)
        self.assertEqual([(e.time, e.kind, e.participants) for e in compiled.events],
                         [(e.time, e.kind, e.participants) for e in reference.events])
        np.testing.assert_array_equal(a.position, b.position)
        np.testing.assert_array_equal(a.velocity, b.velocity)
        self.assertEqual(compiled.ledger.heat_hot.value, reference.ledger.heat_hot.value)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_batched_disc_pair_cases_match_scalar(self):
        cases = ((.2, -1., 1.), (.2, 1000., -1000.), (.1, 0., 0.),
                 (.1, 1., -1.), (0., 0., 0.))
        for distance, va, vb in cases:
            with self.subTest(distance=distance, va=va, vb=vb):
                state = BodyArrays.from_specs([
                    BodySpec((0., 0.), (va, 0.), radius=.05),
                    BodySpec((distance, 0.), (vb, 0.), radius=.05)])
                raw = disc_pairs_toi(state, np.array([[0, 1]], dtype=np.int64),
                                     .1, Tolerances())[0]
                reference = body_pair_toi(state, 0, 1, .1)
                self.assertEqual(int(raw[0]), {"no_collision": 0,
                                               "collision": 1,
                                               "indeterminate": 2}[reference.status.value])
                if reference.contact:
                    self.assertAlmostEqual(raw[1], reference.time, delta=2e-12)
                    np.testing.assert_allclose(raw[3:5], reference.contact.point,
                                               atol=2e-12, rtol=0)
                    np.testing.assert_allclose(raw[5:7], reference.contact.normal,
                                               atol=2e-12, rtol=0)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_batched_disc_pairs_preserve_carnot_chronology(self):
        config = RunConfig(preset="carnot_discs", particles=96, seed=123,
                           pair_search="sweep", numeric_backend="numba")
        batched = load_preset(config)
        scalar = load_preset(config)
        scalar.pair_kernel = "scalar"
        a = batched.advance_to(.05)
        b = scalar.advance_to(.05)
        self.assertEqual([(e.kind, e.participants) for e in batched.events],
                         [(e.kind, e.participants) for e in scalar.events])
        np.testing.assert_allclose([e.time for e in batched.events],
                                   [e.time for e in scalar.events], atol=2e-12, rtol=0)
        np.testing.assert_allclose(a.position, b.position, atol=2e-11, rtol=0)
        np.testing.assert_allclose(a.velocity, b.velocity, atol=2e-11, rtol=0)
        self.assertAlmostEqual(batched.ledger.heat_hot.value,
                               scalar.ledger.heat_hot.value, delta=2e-11)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_controlled_cam_batch_matches_scalar_disc_queries(self):
        sim = load_preset(RunConfig(preset="carnot_discs", particles=24,
                                    numeric_backend="numba"))
        rng = np.random.default_rng(9191)
        wall = sim.world.walls[3]
        state = sim.world.bodies
        for phase in (.2, 1.4, 2.0, 4.6):
            sim.world.mechanism.shaft.phi = phase
            piston_x = wall.point_at(sim.time)[0]
            for a in range(state.n):
                state.pos[a] = (piston_x - state.radius[a] - rng.uniform(.001, .2),
                                rng.uniform(.1, .9))
                state.vel[a] = (rng.uniform(-4, 8), rng.uniform(-3, 3))
            horizon = .03
            speed = wall.speed_bound(sim.time, sim.time + horizon)
            batch = controlled_cam_disc_toi(state, np.ones(state.n, dtype=np.bool_),
                                            sim.world.mechanism, speed, horizon,
                                            sim.tol)
            for a in range(state.n):
                scalar = sim._wall_toi(a, 3, horizon, speed, wall.point_at(sim.time))
                self.assertEqual(int(batch[a, 0]), {"no_collision": 0,
                                                    "collision": 1,
                                                    "indeterminate": 2}[scalar.status.value])
                if scalar.contact:
                    self.assertAlmostEqual(batch[a, 1], scalar.time, delta=2e-11)
                    np.testing.assert_allclose(batch[a, 3:5], scalar.contact.point,
                                               atol=2e-10, rtol=0)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_controlled_cam_batch_preserves_disc_run(self):
        config = RunConfig(preset="carnot_discs", particles=96, seed=125,
                           pair_search="sweep", numeric_backend="numba")
        compiled = load_preset(config)
        scalar = load_preset(config)
        scalar.cam_kernel = "python"
        a = compiled.advance_to(.5)
        b = scalar.advance_to(.5)
        self.assertGreater(sum(bool(e.metadata and e.metadata.get("boundary") == 3)
                               for e in compiled.events), 0)
        self.assertEqual([(e.kind, e.participants) for e in compiled.events],
                         [(e.kind, e.participants) for e in scalar.events])
        np.testing.assert_allclose([e.time for e in compiled.events],
                                   [e.time for e in scalar.events], atol=2e-11, rtol=0)
        np.testing.assert_allclose(a.position, b.position, atol=2e-10, rtol=0)
        np.testing.assert_allclose(a.velocity, b.velocity, atol=2e-9, rtol=0)
        self.assertAlmostEqual(compiled.ledger.heat_hot.value,
                               scalar.ledger.heat_hot.value, delta=2e-9)

    def test_swept_box_search_is_conservative_and_sorted(self):
        rng = np.random.default_rng(919)
        for count in (32, 128):
            pos = rng.uniform(-2, 2, (count, 2))
            vel = rng.normal(0, 6, (count, 2))
            radius = rng.uniform(.002, .07, count)
            for horizon in (0., .03, .2):
                candidates = swept_pairs_sweep(pos, vel, radius, horizon, "python")
                self.assertEqual(candidates, sorted(set(candidates)))
                self.assertTrue(set(candidates) <= set(swept_pairs(pos, vel, radius,
                                                                   horizon)))
                for a in range(count):
                    for b in range(a + 1, count):
                        dp, dv = pos[b] - pos[a], vel[b] - vel[a]
                        speed_sq = float(dv @ dv)
                        closest = (np.clip(-float(dp @ dv) / speed_sq, 0, horizon)
                                   if speed_sq else 0.)
                        if np.linalg.norm(dp + closest * dv) <= radius[a] + radius[b]:
                            self.assertIn((a, b), candidates)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_compiled_swept_box_search_matches_python(self):
        rng = np.random.default_rng(199)
        for count in (32, 200):
            pos = rng.normal(size=(count, 2))
            vel = rng.normal(size=(count, 2))
            radius = rng.uniform(.002, .08, count)
            for horizon in (0., .01, .1):
                self.assertEqual(swept_pairs_sweep(pos, vel, radius, horizon, "python"),
                                 swept_pairs_sweep(pos, vel, radius, horizon, "numba"))

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_swept_box_array_matches_public_pairs_and_solver_order(self):
        rng = np.random.default_rng(279)
        for count in (0, 32, 200):
            pos = rng.normal(size=(count, 2))
            vel = rng.normal(size=(count, 2))
            radius = rng.uniform(.002, .08, count)
            for horizon in (0., .01, .1):
                array = swept_pairs_sweep_array(pos, vel, radius, horizon, "numba")
                self.assertEqual(array.shape[1], 2)
                self.assertEqual([tuple(pair) for pair in array],
                                 swept_pairs_sweep(pos, vel, radius, horizon, "numba"))
                self.assertEqual([tuple(pair) for pair in array],
                                 sorted(set(map(tuple, array))))
        sim = load_preset(RunConfig(preset="carnot_discs", particles=96,
                                    pair_search="sweep", numeric_backend="numba"))
        self.assertIsInstance(sim._candidate_pairs(.01), np.ndarray)
        np.testing.assert_array_equal(
            sim._candidate_pairs(.01),
            swept_pairs_sweep_array(sim.world.bodies.pos, sim.world.bodies.vel,
                                    sim.world.bodies.radius, .01, "numba"))

    def test_sweep_and_grid_have_same_carnot_events(self):
        for preset in ("carnot_discs", "carnot_triangles"):
            config = RunConfig(preset=preset, particles=48, seed=123,
                               max_horizon=.02, numeric_backend="python")
            sweep = load_preset(config)
            sweep.pair_search = "sweep"
            grid = load_preset(config)
            a = sweep.advance_to(.12)
            b = grid.advance_to(.12)
            self.assertEqual([(e.time, e.kind, e.participants) for e in sweep.events],
                             [(e.time, e.kind, e.participants) for e in grid.events])
            np.testing.assert_array_equal(a.position, b.position)
            np.testing.assert_array_equal(a.velocity, b.velocity)
            self.assertEqual(sweep.ledger.work_on.value, grid.ledger.work_on.value)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_compiled_sweep_and_grid_match_large_short_run(self):
        config = RunConfig(preset="carnot_triangles", particles=200, seed=123,
                           numeric_backend="numba")
        sweep = load_preset(config)
        sweep.pair_search = "sweep"
        grid = load_preset(config)
        a = sweep.advance_to(.02)
        b = grid.advance_to(.02)
        self.assertEqual([(e.time, e.kind, e.participants) for e in sweep.events],
                         [(e.time, e.kind, e.participants) for e in grid.events])
        np.testing.assert_array_equal(a.position, b.position)
        np.testing.assert_array_equal(a.velocity, b.velocity)
        self.assertEqual(sweep.ledger.heat_hot.value, grid.ledger.heat_hot.value)
        self.assertEqual(sweep.ledger.work_on.value, grid.ledger.work_on.value)
        self.assertEqual(sweep.max_penetration, grid.max_penetration)

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_batch_polygon_pairs_match_scalar_kernel(self):
        sim = load_preset(RunConfig(preset="carnot_triangles", particles=32,
                                    numeric_backend="numba", seed=173))
        state = sim.world.bodies
        pairs = np.array([(a, b) for a in range(10) for b in range(a + 1, 12)],
                         dtype=np.int64)
        for horizon in (.001, .02, .1):
            batch = polygon_pairs_toi(state, pairs, horizon, sim.tol, 1024)
            for (a, b), raw in zip(pairs, batch):
                scalar = polygon_pair_toi(state, int(a), int(b), horizon,
                                          sim.tol, 1024)
                np.testing.assert_array_equal(raw, np.asarray(scalar, float))

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_batch_pair_kernel_matches_scalar_carnot_chronology(self):
        config = RunConfig(preset="carnot_triangles", particles=96, seed=123,
                           numeric_backend="numba", pair_search="sweep")
        batched = load_preset(config)
        scalar = load_preset(config)
        scalar.pair_kernel = "scalar"
        a = batched.advance_to(.05)
        b = scalar.advance_to(.05)
        self.assertEqual([(e.time, e.kind, e.participants) for e in batched.events],
                         [(e.time, e.kind, e.participants) for e in scalar.events])
        np.testing.assert_array_equal(a.position, b.position)
        np.testing.assert_array_equal(a.velocity, b.velocity)
        self.assertEqual(batched.ledger.heat_hot.value, scalar.ledger.heat_hot.value)
        self.assertEqual(batched.ledger.work_on.value, scalar.ledger.work_on.value)

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

    @unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
    def test_bounded_moving_wall_queries_through_two_complete_cycles(self):
        end = 4*math.pi/1.5
        for seed in (123,124):
            with self.subTest(seed=seed):
                config=RunConfig(preset="carnot_triangles",particles=32,
                                 seed=seed,shaft_speed=1.5,max_horizon=.05,
                                 pair_search="sweep",numeric_backend="numba")
                bounded=load_preset(config)
                checkpoint=None
                for cycle in range(2):
                    for branch in range(4):
                        midpoint=(cycle+(branch+.5)/4)*2*math.pi/1.5
                        bounded.advance_to(midpoint)
                        if seed==123 and cycle==1 and branch==0:
                            checkpoint=bounded.checkpoint()
                        local=bounded._earliest(.05)
                        bounded.wall_search="all"
                        unpruned=bounded._earliest(.05)
                        bounded.wall_search="bounded"
                        self.assertEqual(local.status,unpruned.status)
                        self.assertEqual(local.reason,unpruned.reason)
                        if local.time is not None:
                            self.assertAlmostEqual(local.time,unpruned.time,
                                                   delta=2e-11)
                        self.assertEqual(len(local.contacts),len(unpruned.contacts))
                        for a,b in zip(local.contacts,unpruned.contacts):
                            self.assertEqual((a.a,a.b,a.boundary,a.feature_a,a.feature_b),
                                             (b.a,b.b,b.boundary,b.feature_a,b.feature_b))
                            np.testing.assert_allclose(a.point,b.point,
                                                       atol=2e-10,rtol=0)
                            np.testing.assert_allclose(a.normal,b.normal,
                                                       atol=2e-10,rtol=0)
                bounded.advance_to(end)

                self.assertEqual(len(bounded.cycle_markers),3)
                self.assertEqual(sum(event.kind=="branch_transition"
                                     for event in bounded.events),8)
                self.assertGreater(sum(event.metadata is not None and
                                       event.metadata.get("boundary")==3
                                       for event in bounded.events),0)
                self.assertEqual(bounded.ccd_failures,0)
                self.assertLessEqual(bounded.max_penetration,bounded.tol.geometry)
                self.assertLess(abs(bounded.snapshot().energy_residual),1e-10)
                if seed==123:
                    direct=load_preset(config)
                    direct.advance_to(end)
                    expected_events=[event.as_dict() for event in direct.events]
                    self.assertEqual([event.as_dict() for event in bounded.events],
                                     expected_events)
                    np.testing.assert_array_equal(bounded.world.bodies.pos,
                                                  direct.world.bodies.pos)
                    self.assertEqual(bounded.ledger,direct.ledger)
                    bounded.restore(checkpoint)
                    bounded.advance_to(end)
                    self.assertEqual([event.as_dict() for event in bounded.events],
                                     expected_events)
                    np.testing.assert_array_equal(bounded.world.bodies.pos,
                                                  direct.world.bodies.pos)
                    self.assertEqual(bounded.ledger,direct.ledger)

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
    def test_batched_triangle_wall_matches_scalar_edge_contacts(self):
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
                                      numeric_backend="numba", wall_search="bounded")
                scalar = Simulation(World(BodyArrays.from_specs(specs), [wall]),
                                    numeric_backend="numba", wall_search="bounded",
                                    wall_kernel="python")
                a = compiled._earliest(horizon)
                b = scalar._earliest(horizon)
                self.assertEqual(a.status, b.status)
                self.assertEqual(a.reason, b.reason)
                if a.contact is not None:
                    self.assertAlmostEqual(a.time, b.time, delta=2e-11)
                    np.testing.assert_allclose(a.contact.point, b.contact.point,
                                               atol=2e-10, rtol=0)
                    self.assertEqual(a.contact.feature_a, b.contact.feature_a)

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
