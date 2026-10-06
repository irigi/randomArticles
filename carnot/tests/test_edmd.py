"""Compiled event-driven disc kernel (osmosis plan, milestone 1)."""

import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.core import edmd
from microthermo.core.boundaries import Wall
from microthermo.core.state import BodyArrays, BodySpec, Shape
from microthermo.experiments import build_preset
from microthermo.experiments.common import box_walls
from microthermo.measurements.live import LiveInstruments
from microthermo.runner.edmd_simulation import EdmdSimulation, edmd_unsupported_reason
from microthermo.runner.simulation import NumericalFailure, World


def world_of(specs, width=2., height=1.):
    return World(BodyArrays.from_specs(specs), box_walls(width, height))


def gas_box(particles=32, seed=123, engine="edmd"):
    return load_preset(RunConfig(preset="gas_box", particles=particles, seed=seed,
                                 engine=engine))


def gas_box_world(particles):
    return build_preset(RunConfig(preset="gas_box", particles=particles),
                        np.random.default_rng(123))


@unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
class EdmdAnalyticTests(unittest.TestCase):
    def test_head_on_pair_time_and_exchange(self):
        sim = EdmdSimulation(world_of([BodySpec((0.5, .5), (1., 0.), radius=.05),
                                       BodySpec((1.5, .5), (-1., 0.), radius=.05)]))
        sim.step_collision()
        self.assertAlmostEqual(sim.time, .45, delta=1e-15)
        event = sim.events[-1]
        self.assertEqual((event.kind, event.participants), ("elastic", (0, 1)))
        np.testing.assert_allclose(sim.world.bodies.vel, [[-1., 0.], [1., 0.]], atol=1e-15)
        np.testing.assert_allclose(event.impulse, (2., 0.))

    def test_unequal_masses_glancing_conserves_energy_and_momentum(self):
        sim = EdmdSimulation(world_of([BodySpec((0.5, .5), (1., 0.), mass=1., radius=.05),
                                       BodySpec((1.0, .53), (0., 0.), mass=3., radius=.05)]))
        p0, e0 = sim.world.bodies.momentum().copy(), sim.energy()
        sim.step_collision()
        self.assertEqual(sim.events[-1].kind, "elastic")
        np.testing.assert_allclose(sim.world.bodies.momentum(), p0, atol=1e-15)
        self.assertAlmostEqual(sim.energy(), e0, delta=1e-15*e0)
        d = np.linalg.norm(sim.world.bodies.pos[1] - sim.world.bodies.pos[0])
        self.assertAlmostEqual(d, .1, delta=1e-14)

    def test_specular_wall_hit_time_and_reflection(self):
        sim = EdmdSimulation(world_of([BodySpec((1., .5), (2., 1.), radius=.05)]))
        sim.step_collision()
        self.assertAlmostEqual(sim.time, .45, delta=1e-15)   # top wall, gap .45
        self.assertEqual(sim.events[-1].kind, "specular")
        np.testing.assert_allclose(sim.world.bodies.vel[0], [2., -1.])

    def test_step_collision_commits_exactly_one_event(self):
        sim = gas_box(16)
        for k in range(1, 30):
            sim.step_collision()
            self.assertEqual(sim.event_count, k)
            self.assertEqual(len(sim.events), k)

    def test_overlap_is_a_numerical_failure_with_checkpoint(self):
        sim = EdmdSimulation(world_of([BodySpec((0.5, .5), (1., 0.), radius=.05),
                                       BodySpec((1.5, .5), (-1., 0.), radius=.05)]))
        sim.pos[1] = (0.56, .5)        # corrupt the state: overlapping discs
        with self.assertRaises(NumericalFailure):
            sim.advance_to(.1)
        with tempfile.TemporaryDirectory() as tmp:
            path = sim.save_failure(Path(tmp)/"failure")
            self.assertTrue((path/"diagnostic.json").exists())
            self.assertTrue((path/"state.npz").exists())

    def test_unsupported_worlds_are_rejected(self):
        triangle = World(BodyArrays.from_specs([BodySpec((1., .5), (0., 0.), radius=.05,
                                                         shape=Shape.TRIANGLE)]),
                         box_walls(2., 1.))
        self.assertIn("discs", edmd_unsupported_reason(triangle))
        moving = box_walls(2., 1.)
        moving[1] = Wall(np.array([2., .5]), np.array([-1., 0.]), 1.,
                         velocity=np.array([-.1, 0.]), name="piston")
        disc = BodyArrays.from_specs([BodySpec((1., .5), (0., 0.), radius=.05)])
        self.assertIn("moves", edmd_unsupported_reason(World(disc, moving)))
        with self.assertRaises(ValueError):
            load_preset(RunConfig(preset="carnot_discs", engine="edmd"))
        with self.assertRaises(ValueError):
            RunConfig(engine="warp").validate()

    def test_engine_default_keeps_old_config_digests(self):
        self.assertEqual(RunConfig().digest(), RunConfig(engine="reference").digest())
        self.assertNotEqual(RunConfig().digest(), RunConfig(engine="edmd").digest())


@unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
class EdmdDynamicsTests(unittest.TestCase):
    def test_closed_specular_box_conserves_energy_and_never_overlaps(self):
        rng = np.random.default_rng(7)
        specs = []
        while len(specs) < 300:
            p = rng.uniform(.03, [3.97, 1.97])
            if all(np.hypot(*(p - q.position)) > .062 for q in specs):
                specs.append(BodySpec(tuple(p), tuple(rng.normal(0, 1, 2)), radius=.03))
        sim = EdmdSimulation(world_of(specs, 4., 2.))
        e0 = sim.energy()
        for k in range(1, 21):
            sim.advance_to(k*.5)          # each sample checks the minimum gap
        self.assertGreater(sim.event_count, 10000)
        self.assertLessEqual(abs(sim.energy() - e0), 1e-12*e0)
        self.assertLessEqual(abs(sim.snapshot().energy_residual), 1e-12*e0)
        self.assertLessEqual(sim.max_penetration, 1e-12)

    def test_thermal_wall_first_law_and_flux_distribution(self):
        sim = gas_box(64)
        sim.advance_to(100.)
        snap = sim.snapshot()
        self.assertLessEqual(abs(snap.energy_residual), 1e-11*snap.energy)
        out = np.array([e.metadata["outgoing_normal_mode_energy"] for e in sim.events
                        if e.kind == "hot"])
        # Flux-weighted refresh: the outgoing normal-mode energy is Exp(T).
        self.assertGreater(len(out), 1000)
        self.assertAlmostEqual(np.mean(out), 1., delta=4/math.sqrt(len(out)))
        self.assertAlmostEqual(np.mean(out**2), 2., delta=25/math.sqrt(len(out)))
        heat = sum(e.heat_into_system for e in sim.events)
        self.assertAlmostEqual(heat, sim.ledger.heat_hot.value, delta=1e-10)

    def test_gas_relaxes_to_the_wall_temperature(self):
        sim = gas_box(64)
        sim._bring_to_now()
        sim.vel *= math.sqrt(2.)           # start hot; one thermal wall at T = 1
        sim._changed_velocities()
        temperatures = [sim.advance_to(100. + k).translational_temperature
                        for k in range(1, 201)]
        self.assertAlmostEqual(np.mean(temperatures), 1., delta=.06)

    def test_reversal_retraces_a_specular_run(self):
        rng = np.random.default_rng(3)
        specs = [BodySpec((x, y), tuple(rng.normal(0, 1, 2)), radius=.04)
                 for x, y in [(.3, .3), (.9, .4), (1.5, .6), (.6, .75), (1.2, .2), (1.7, .3)]]
        sim = EdmdSimulation(world_of(specs))
        start = sim.world.bodies.pos.copy()
        sim.advance_to(1.5)
        self.assertGreater(sim.event_count, 10)
        sim.apply_command({"name": "reverse_particle_velocities"})
        sim.advance_to(3.)
        np.testing.assert_allclose(sim.world.bodies.pos, start, atol=1e-8)


@unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
class EdmdReproducibilityTests(unittest.TestCase):
    def test_sample_cadence_does_not_change_the_trajectory(self):
        coarse, fine = gas_box(48), gas_box(48)
        coarse.advance_to(6.)
        for k in range(1, 601):
            fine.advance_to(k*.01)
        self.assertEqual(coarse.event_count, fine.event_count)
        np.testing.assert_array_equal(coarse.world.bodies.pos, fine.world.bodies.pos)
        np.testing.assert_array_equal(coarse.world.bodies.vel, fine.world.bodies.vel)

    def test_checkpoint_restore_is_bitwise_reproducible(self):
        sim = gas_box(48)
        sim.advance_to(2.)
        cp = sim.checkpoint()
        first = sim.advance_to(5.)
        events = [(e.time, e.kind, e.participants) for e in sim.events]
        sim.restore(cp)
        self.assertEqual(sim.time, 2.)
        second = sim.advance_to(5.)
        np.testing.assert_array_equal(first.position, second.position)
        np.testing.assert_array_equal(first.velocity, second.velocity)
        self.assertEqual([(e.time, e.kind, e.participants) for e in sim.events], events)

    def test_buffer_and_log_sizes_do_not_change_the_run(self):
        a = gas_box(48)
        b = EdmdSimulation(gas_box_world(48), 123, random_block=7, log_capacity=5)
        a.advance_to(5.)
        b.advance_to(5.)
        self.assertEqual(a.event_count, b.event_count)
        np.testing.assert_array_equal(a.world.bodies.pos, b.world.bodies.pos)
        self.assertEqual([e.time for e in a.events], [e.time for e in b.events])

    def test_matches_reference_engine_until_chaos_amplifies_wall_tolerance(self):
        ref, new = gas_box(8, engine="reference"), gas_box(8)
        ref.advance_to(4.)
        new.advance_to(4.)
        pairs = list(zip(ref.events, new.events))[:40]
        self.assertEqual(len(pairs), 40)
        # Disc-disc times are exact quadratics in both engines; reference wall
        # contacts are found within a 1e-10 gap, which chaos then amplifies.
        for a, b in pairs:
            self.assertEqual((a.kind, tuple(a.participants)), (b.kind, tuple(b.participants)))
            self.assertEqual((a.metadata or {}).get("boundary"),
                             (b.metadata or {}).get("boundary"))
            self.assertAlmostEqual(b.time, a.time, delta=1e-6)
        # Identical uniform draws give the same flux-weighted outgoing energies.
        hot = [(a, b) for a, b in pairs if a.kind == "hot"]
        self.assertTrue(hot)
        for a, b in hot:
            expected = a.metadata["outgoing_normal_mode_energy"]
            self.assertAlmostEqual(b.metadata["outgoing_normal_mode_energy"], expected,
                                   delta=1e-12*expected)

    def test_grid_and_single_cell_agree_until_chaos(self):
        big = gas_box(64)
        one = EdmdSimulation(gas_box_world(64), 123)
        one.ip[:] = 1
        one.fp[edmd.FP_CX], one.fp[edmd.FP_CY] = 2., 1.
        one.head = np.full(1, -1, np.int64)
        one.cps, one.cpi = np.zeros(2, np.int64), np.zeros(0, np.int64)
        one.reg = np.zeros((1, 0), np.int8)
        edmd.build_grid(one._bodies(), one._grid(), one._hosts(), 0.)
        one._rebuild()
        self.assertGreater(big.ip[0]*big.ip[1], 100)
        big.advance_to(.5)
        one.advance_to(.5)
        pairs = list(zip(big.events, one.events))[:100]
        self.assertEqual(len(pairs), 100)
        for a, b in pairs:
            self.assertEqual((a.kind, a.participants), (b.kind, b.participants))
            self.assertAlmostEqual(b.time, a.time, delta=1e-9)

    def test_live_instruments_observe_the_edmd_engine(self):
        sim = gas_box(32)
        instruments = LiveInstruments()
        for k in range(1, 6):
            instruments.observe(sim, sim.advance_to(k*.2))
        self.assertGreater(sim.event_count, 0)


if __name__ == "__main__":
    unittest.main()
