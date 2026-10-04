"""The benchmark matrix records independently measured, reusable runs."""

import json
from pathlib import Path
import tempfile
import unittest

from microthermo.measurements.performance_matrix import main


class PerformanceMatrixTests(unittest.TestCase):
    def test_two_counts_record_metadata_and_throughput(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "nested" / "matrix.json"
            result = main(["--presets", "carnot_discs", "--counts", "16,32",
                           "--duration", ".01", "--numeric-backend", "python",
                           "--output", str(output)])
            self.assertEqual(result, 0)
            report = json.loads(output.read_text())
            self.assertTrue(report["complete"])
            self.assertEqual(report["schema_version"], 1)
            self.assertEqual([row["particles"] for row in report["runs"]], [16, 32])
            for row in report["runs"]:
                self.assertEqual(row["numeric_backend"], "python")
                self.assertEqual(row["seed"], 123)
                self.assertGreater(row["process_seconds"], row["solver_seconds"])
                self.assertGreater(row["peak_process_rss_kib"], 0)
                self.assertAlmostEqual(row["wall_seconds_per_physical_second"],
                                       row["solver_seconds"] / row["physical_seconds"])


if __name__ == "__main__":
    unittest.main()
