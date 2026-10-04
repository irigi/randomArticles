"""Small end-to-end check for the long-run measurement command."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from microthermo.config import RunConfig
from microthermo.core.numeric import numba_available
from microthermo.measurements import triangle_cycle_health
from microthermo.measurements.triangle_cycle_health import measure


@unittest.skipUnless(numba_available(), "Numba optional dependency not installed")
class TriangleCycleHealthTests(unittest.TestCase):
    def test_midpoint_queries_preserve_direct_replay(self):
        config = RunConfig(preset="carnot_triangles", particles=12, seed=123,
                           shaft_speed=1.5, pair_search="sweep",
                           numeric_backend="numba")
        result = measure(config, branches=4, cadence=.083,
                         same_state_queries=True)
        self.assertEqual(len(result["same_state_queries"]), 4)
        self.assertTrue(all(query["matches"]
                            for query in result["same_state_queries"]))
        self.assertTrue(result["sampling"]["matches_direct"])
        self.assertTrue(result["sampling"]["replay_matches_direct"])
        self.assertEqual(result["ccd_failures"], 0)

    def test_branch_checkpoint_resumes_exact_history(self):
        config = RunConfig(preset="carnot_triangles", particles=12, seed=124,
                           shaft_speed=1.5, pair_search="sweep",
                           numeric_backend="numba")
        with TemporaryDirectory() as directory:
            checkpoint = Path(directory) / "branch.pkl.gz"
            writer = triangle_cycle_health.save_branch_checkpoint

            def interrupt_after_first(path, payload):
                writer(path, payload)
                if payload["completed_branch"] == 1:
                    raise RuntimeError("simulated interruption")

            with patch.object(triangle_cycle_health, "save_branch_checkpoint",
                              side_effect=interrupt_after_first):
                with self.assertRaisesRegex(RuntimeError, "simulated interruption"):
                    measure(config, branches=2, cadence=None,
                            same_state_queries=True,
                            branch_checkpoint=checkpoint)
            resumed = measure(config, branches=2, cadence=None,
                              same_state_queries=True,
                              branch_checkpoint=checkpoint)
            direct = measure(config, branches=2, cadence=None,
                             same_state_queries=True)
            for name in ("events", "event_kinds", "completed_cycles",
                         "ccd_refinements", "ccd_failures", "max_penetration",
                         "energy_residual", "ledger", "branches",
                         "same_state_queries"):
                self.assertEqual(resumed[name], direct[name], name)


if __name__ == "__main__":
    unittest.main()
