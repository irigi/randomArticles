"""Osmosis replays (archive version 3) and MP4 export (osmosis plan, milestone 7)."""

import dataclasses
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import numpy as np

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
try:
    from PySide6 import QtWidgets
except ImportError:
    QtWidgets = None

from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.core import edmd
from microthermo.io.replay import ReplayReader, precalculate_replay
from microthermo.measurements.osmosis import OsmosisInstruments

FPS = 4.
CONFIG = RunConfig(preset="osmosis", particles=60, hosts=4, duration=50., seed=5,
                   temperature_schedule=((25., 2.),))


def same(a, b):
    """Dataclass equality that treats NaN (an empty well's mu) as equal to NaN."""
    def fields(x):
        return tuple("nan" if isinstance(v, float) and math.isnan(v) else v
                     for v in dataclasses.astuple(x))
    return fields(a) == fields(b)


@unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
class OsmosisReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.TemporaryDirectory()
        cls.path = precalculate_replay(CONFIG, Path(cls.directory.name)/"replay", fps=FPS,
                                       chunk_frames=16)

    @classmethod
    def tearDownClass(cls):
        cls.directory.cleanup()

    def test_recorded_osmosis_frames_match_the_live_instruments(self):
        reader = ReplayReader(self.path)
        self.assertEqual(reader.manifest["version"], 3)
        self.assertTrue(reader.has_osmosis)
        sim = load_preset(CONFIG)
        live = OsmosisInstruments()
        checked, since = 0, set()
        for index in range(len(reader)):
            snapshot = sim.advance_to(index/FPS) if index else sim.snapshot()
            expected = live.observe(sim, snapshot)
            if index % 7 and index != len(reader) - 1:
                continue
            recorded = reader.osmosis_frame(index)
            np.testing.assert_array_equal(reader.frame(index).membership, sim.memberships)
            self.assertTrue(same(recorded.current, expected.current))
            self.assertEqual(len(recorded.history), len(expected.history))
            self.assertTrue(all(map(same, recorded.history, expected.history)))
            self.assertEqual(recorded.reference, expected.reference)
            self.assertEqual(recorded.averages.keys(), expected.averages.keys())
            for key, value in expected.averages.items():
                self.assertAlmostEqual(recorded.averages[key], value, delta=1e-12*abs(value))
            np.testing.assert_array_equal(recorded.occupancy_histogram,
                                          expected.occupancy_histogram)
            checked += 1
            since.add(recorded.averages.get("since"))
        self.assertGreater(checked, 20)
        # Averaging starts at t = 20 and again 20 s after the step at t = 25.
        self.assertEqual(since, {None, 20., 45.})
        self.assertIn("osmotic_pressure", recorded.averages)
        self.assertEqual(recorded.current.temperature, 2.)

    def test_membership_between_frames_follows_step_events(self):
        reader = ReplayReader(self.path)
        reader.prepare_playback()
        try:
            changed = 0
            for index in range(len(reader) - 1):
                left, right = reader.frame(index), reader.frame(index + 1)
                if np.array_equal(left.membership, right.membership):
                    continue
                changed += 1
                middle, _ = reader.sample(.5*(left.time + right.time))
                late, _ = reader.sample(right.time - 1e-9)
                # Every disc shows its left or its right membership.
                self.assertTrue(np.all((middle.membership == left.membership) |
                                       (middle.membership == right.membership)))
                np.testing.assert_array_equal(late.membership, right.membership)
            self.assertGreater(changed, 5)
        finally:
            reader.close()

    @unittest.skipIf(QtWidgets is None, "PySide6 is not installed")
    def test_replay_window_shows_the_osmosis_panel(self):
        from microthermo.ui.replay_window import ReplayWindow

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        window = ReplayWindow(self.path)
        try:
            window.show()
            app.processEvents()
            window.slider.setValue(9000)
            app.processEvents()
            panel = window.lab.osmosis
            self.assertTrue(panel.isVisibleTo(window))
            self.assertFalse(window.lab.diagrams.isVisibleTo(window))
            self.assertIn("osmotic pressure", panel.summary.text())
            self.assertGreater(len(panel.left_line.getData()[0]), 10)
            self.assertEqual([t for t, _ in panel._steps], [25.])
            self.assertIn("wall T=2.000", window.status_line.text())
            self.assertIsNotNone(window.view.snapshot.membership)
        finally:
            window.close()

    @unittest.skipIf(QtWidgets is None or shutil.which("ffmpeg") is None,
                     "PySide6 or ffmpeg is not available")
    def test_video_export_writes_every_frame(self):
        from microthermo.ui.video import export_video

        output = Path(self.directory.name)/"replay.mp4"
        calls = []
        export_video(self.path, output, fps=10., speed=5., start=20., end=26.,
                     size=(480, 300), progress=lambda k, n: calls.append((k, n)))
        self.assertEqual(calls[-1], (13, 13))
        probe = subprocess.run(
            ["ffprobe", "-v", "error", "-count_frames", "-show_entries",
             "stream=width,height,nb_read_frames", "-of", "csv=p=0", str(output)],
            capture_output=True, text=True, check=True)
        self.assertEqual(probe.stdout.strip(), "480,300,13")
        with self.assertRaises(ValueError):
            export_video(self.path, output, size=(481, 300))


if __name__ == "__main__":
    unittest.main()
