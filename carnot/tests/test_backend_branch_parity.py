"""Full-branch checks for the optional compiled triangle backend."""

import math
import unittest

import numpy as np

from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.core.numeric import numba_available


COUNTERS = ("heat_hot", "heat_cold", "heat_other", "work_on",
            "load_output", "piston_work_on_gas", "support_impulse_x",
            "support_impulse_y")


@unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
class BackendBranchParityTests(unittest.TestCase):
    def make_sim(self, seed, backend, particles=12, shaft_speed=1.5):
        config = RunConfig(preset="carnot_triangles", particles=particles, seed=seed,
                           shaft_speed=shaft_speed, max_horizon=.05,
                           pair_search="sweep", numeric_backend=backend)
        return load_preset(config)

    def assert_parity(self, compiled, reference):
        self.assertEqual(len(compiled.events), len(reference.events))
        for a, b in zip(compiled.events, reference.events):
            self.assertEqual((a.kind, a.participants), (b.kind, b.participants))
            self.assertAlmostEqual(a.time, b.time, delta=5e-11)
            np.testing.assert_allclose(a.impulse, b.impulse, atol=4e-7, rtol=0)
            for field in ("energy_before", "energy_after", "heat_into_system",
                          "work_on_system"):
                self.assertAlmostEqual(getattr(a, field), getattr(b, field),
                                       delta=4e-7)
            self.assertEqual(set(a.metadata or ()), set(b.metadata or ()))
            for key in a.metadata or ():
                left, right = a.metadata[key], b.metadata[key]
                if isinstance(left, (int, float)) and not isinstance(left, bool):
                    self.assertAlmostEqual(left, right, delta=4e-7)
                else:
                    self.assertEqual(left, right)
        for field, tol in (("pos", 1e-8), ("vel", 4e-7),
                           ("angle", 1e-6), ("omega", 2e-5)):
            np.testing.assert_allclose(getattr(compiled.world.bodies, field),
                                       getattr(reference.world.bodies, field),
                                       atol=tol, rtol=0)
        for field in COUNTERS:
            self.assertAlmostEqual(getattr(compiled.ledger, field).value,
                                   getattr(reference.ledger, field).value,
                                   delta=2e-8)
        self.assertEqual(compiled.ccd_failures, reference.ccd_failures)
        self.assertEqual(compiled.failure_diagnostic, reference.failure_diagnostic)
        self.assertLessEqual(compiled.max_penetration, compiled.tol.geometry)
        self.assertLessEqual(reference.max_penetration, reference.tol.geometry)
        self.assertAlmostEqual(compiled.snapshot().energy_residual,
                               reference.snapshot().energy_residual, delta=2e-8)

    def test_complete_hot_branch_across_seeds_and_sampling_cadences(self):
        end = math.pi / 3  # one quarter turn at prescribed speed 1.5
        for seed in (123, 124):
            with self.subTest(seed=seed):
                direct = self.make_sim(seed, "numba")
                reference = self.make_sim(seed, "python")
                direct.advance_to(end)
                reference.advance_to(end)
                self.assertEqual(direct.world.mechanism.branch,
                                 "adiabatic_expansion")
                self.assertEqual(sum(e.kind == "branch_transition"
                                     for e in direct.events), 1)
                self.assert_parity(direct, reference)

                for backend in ("numba", "python"):
                    sampled = self.make_sim(seed, backend)
                    # Two cadences cross different event intervals. Restoring
                    # mid-branch must reproduce the same end state and ledger.
                    cadence = .083 if seed == 123 else .137
                    checkpoint = None
                    for target in np.arange(cadence, end, cadence):
                        sampled.advance_to(float(target))
                        if checkpoint is None and target > end / 3:
                            checkpoint = sampled.checkpoint()
                    sampled.advance_to(end)
                    baseline = direct if backend == "numba" else reference
                    self.assertEqual(sampled.events, baseline.events)
                    np.testing.assert_array_equal(sampled.world.bodies.pos,
                                                  baseline.world.bodies.pos)
                    self.assertEqual(sampled.ledger, baseline.ledger)
                    sampled.restore(checkpoint)
                    sampled.advance_to(end)
                    self.assertEqual(sampled.events, baseline.events)
                    np.testing.assert_array_equal(sampled.world.bodies.pos,
                                                  baseline.world.bodies.pos)
                    self.assertEqual(sampled.ledger, baseline.ledger)

    def test_same_state_queries_match_after_full_cycle_trajectories_drift(self):
        # Event times in independent trajectories separate after the first
        # branch. Re-query the observed sensitive contacts from identical
        # checkpoints to distinguish trajectory amplification from a different
        # collision answer for the same input state.
        for time, horizon in ((1.24, .05), (1.75, .05), (1.96, .1)):
            with self.subTest(time=time):
                compiled = self.make_sim(123, "numba")
                compiled.advance_to(time)
                reference = self.make_sim(123, "python")
                reference.restore(compiled.checkpoint())
                # Both solvers must construct a new query from this checkpoint.
                compiled._pending = None
                reference._pending = None
                a = compiled._earliest(horizon)
                b = reference._earliest(horizon)
                self.assertEqual((a.status, a.reason), (b.status, b.reason))
                self.assertAlmostEqual(a.time, b.time, delta=2e-11)
                self.assertEqual((a.contact.a, a.contact.b, a.contact.boundary,
                                  a.contact.feature_a, a.contact.feature_b),
                                 (b.contact.a, b.contact.b, b.contact.boundary,
                                  b.contact.feature_a, b.contact.feature_b))
                np.testing.assert_allclose(a.contact.point, b.contact.point,
                                           atol=2e-10, rtol=0)
                np.testing.assert_allclose(a.contact.normal, b.contact.normal,
                                           atol=2e-10, rtol=0)

    def test_full_cycle_local_queries_and_sampling_health_at_32_triangles(self):
        duration = 2 * math.pi / 1.5
        branch_names = ("hot", "adiabatic_expansion", "cold",
                        "adiabatic_compression")
        for seed, cadence in ((123, .083), (124, .137)):
            with self.subTest(seed=seed):
                # Increase density of contacts while preserving the configured
                # dilute occupied-area fraction.
                direct = self.make_sim(seed, "numba", particles=32)
                reference = self.make_sim(seed, "python", particles=32)
                for branch_index, branch_name in enumerate(branch_names):
                    time = (branch_index + .5) * math.pi / (2 * 1.5)
                    direct.advance_to(time)
                    self.assertEqual(direct.world.mechanism.branch, branch_name)
                    reference.restore(direct.checkpoint())
                    compiled_query = direct._earliest(.05)
                    python_query = reference._earliest(.05)
                    self.assertEqual(compiled_query.status, python_query.status)
                    self.assertEqual(compiled_query.reason, python_query.reason)
                    self.assertIsNotNone(compiled_query.contact)
                    self.assertAlmostEqual(compiled_query.time,
                                           python_query.time, delta=2e-11)
                    self.assertEqual(len(compiled_query.contacts),
                                     len(python_query.contacts))
                    for a, b in zip(compiled_query.contacts,
                                    python_query.contacts):
                        self.assertEqual((a.a, a.b, a.boundary,
                                          a.feature_a, a.feature_b),
                                         (b.a, b.b, b.boundary,
                                          b.feature_a, b.feature_b))
                        np.testing.assert_allclose(a.point, b.point,
                                                   atol=2e-10, rtol=0)
                        np.testing.assert_allclose(a.normal, b.normal,
                                                   atol=2e-10, rtol=0)
                direct.advance_to(duration)
                sampled = self.make_sim(seed, "numba", particles=32)
                checkpoint = None
                for target in np.arange(cadence, duration, cadence):
                    sampled.advance_to(float(target))
                    if checkpoint is None and target > duration / 2:
                        checkpoint = sampled.checkpoint()
                sampled.advance_to(duration)
                self.assertEqual(sampled.events, direct.events)
                np.testing.assert_array_equal(sampled.world.bodies.pos,
                                              direct.world.bodies.pos)
                np.testing.assert_array_equal(sampled.world.bodies.vel,
                                              direct.world.bodies.vel)
                self.assertEqual(sampled.ledger, direct.ledger)
                self.assertEqual(sampled.cycle_markers, direct.cycle_markers)
                self.assertEqual(len(direct.cycle_markers), 2)
                self.assertEqual(direct.ccd_failures, 0)
                self.assertLessEqual(direct.max_penetration, direct.tol.geometry)
                self.assertLess(abs(direct.snapshot().energy_residual), 1e-9)
                sampled.restore(checkpoint)
                sampled.advance_to(duration)
                self.assertEqual(sampled.events, direct.events)
                np.testing.assert_array_equal(sampled.world.bodies.pos,
                                              direct.world.bodies.pos)
                self.assertEqual(sampled.ledger, direct.ledger)

    def test_default_speed_cycle_health_and_local_parity(self):
        speed = .15
        duration = 2 * math.pi / speed
        branches = ("hot", "adiabatic_expansion", "cold",
                    "adiabatic_compression")
        for seed, cadence in ((123, .83), (124, 1.37), (125, .71)):
            with self.subTest(seed=seed):
                direct = self.make_sim(seed, "numba", 32, speed)
                reference = self.make_sim(seed, "python", 32, speed)
                for index, branch in enumerate(branches):
                    time = (index + .5) * math.pi / (2 * speed)
                    direct.advance_to(time)
                    self.assertEqual(direct.world.mechanism.branch, branch)
                    reference.restore(direct.checkpoint())
                    a, b = direct._earliest(.05), reference._earliest(.05)
                    self.assertEqual((a.status, a.reason), (b.status, b.reason))
                    self.assertEqual(len(a.contacts), len(b.contacts))
                    if a.contact is not None:
                        self.assertAlmostEqual(a.time, b.time, delta=2e-11)
                    for left, right in zip(a.contacts, b.contacts):
                        self.assertEqual((left.a, left.b, left.boundary,
                                          left.feature_a, left.feature_b),
                                         (right.a, right.b, right.boundary,
                                          right.feature_a, right.feature_b))
                        np.testing.assert_allclose(left.point, right.point,
                                                   atol=2e-10, rtol=0)
                        np.testing.assert_allclose(left.normal, right.normal,
                                                   atol=2e-10, rtol=0)
                direct.advance_to(duration)
                sampled = self.make_sim(seed, "numba", 32, speed)
                checkpoint = None
                for target in np.arange(cadence, duration, cadence):
                    sampled.advance_to(float(target))
                    if checkpoint is None and target > duration / 2:
                        checkpoint = sampled.checkpoint()
                sampled.advance_to(duration)
                self.assertGreater(len(direct.events), 4000)
                self.assertEqual(len(direct.cycle_markers), 2)
                self.assertEqual(direct.ccd_failures, 0)
                self.assertIsNone(direct.failure_diagnostic)
                self.assertLessEqual(direct.max_penetration,
                                     direct.tol.geometry)
                self.assertLess(abs(direct.snapshot().energy_residual), 1e-9)
                self.assertEqual(sampled.events, direct.events)
                np.testing.assert_array_equal(sampled.world.bodies.pos,
                                              direct.world.bodies.pos)
                np.testing.assert_array_equal(sampled.world.bodies.vel,
                                              direct.world.bodies.vel)
                self.assertEqual(sampled.ledger, direct.ledger)
                self.assertEqual(sampled.cycle_markers, direct.cycle_markers)
                sampled.restore(checkpoint)
                sampled.advance_to(duration)
                self.assertEqual(sampled.events, direct.events)
                np.testing.assert_array_equal(sampled.world.bodies.pos,
                                              direct.world.bodies.pos)
                self.assertEqual(sampled.ledger, direct.ledger)
