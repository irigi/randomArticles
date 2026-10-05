"""Offscreen smoke test for solver-free replay controls."""

import gc
import os
from pathlib import Path
import tempfile
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
try:
    from PySide6 import QtWidgets
except ImportError:
    QtWidgets = None

from microthermo.config import RunConfig
from microthermo.io.replay import precalculate_replay


@unittest.skipIf(QtWidgets is None,"PySide6 is not installed")
class ReplayUiTests(unittest.TestCase):
    def test_seek_and_playback_use_recorded_frames(self):
        from microthermo.ui.replay_window import ReplayWindow

        app=QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        with tempfile.TemporaryDirectory() as directory:
            config=RunConfig(preset="carnot_discs",particles=4,seed=123,
                             duration=.2)
            path=precalculate_replay(config,Path(directory)/"replay",fps=20.,
                                     chunk_frames=2)
            window=ReplayWindow(path)
            try:
                window.show()
                app.processEvents()
                window.slider.setValue(5000)
                self.assertAlmostEqual(window.view.snapshot.time,.1)
                self.assertIn("recorded statistics at t=0.100",window.status_line.text())
                self.assertTrue(window.lab.diagrams.isVisibleTo(window))
                self.assertNotEqual(window.lab.readout.values["Area"].text(),"—")
                self.assertGreater(len(window.lab.readout.temp_trans.getData()[0]),1)
                window.speed.setCurrentIndex(4)
                self.assertEqual(window.playback_speed,4.)
                # Play slowly so the short archive cannot end during the
                # first (slower) paint of the plots.
                window.speed.setCurrentIndex(0)
                gc_before=gc.isenabled()
                window.play.setChecked(True)
                app.processEvents()
                self.assertTrue(window.timer.isActive())
                self.assertFalse(gc.isenabled())
                window.play.setChecked(False)
                self.assertFalse(window.timer.isActive())
                self.assertEqual(gc.isenabled(),gc_before)
                self.assertGreater(window.grab().toImage().width(),0)
            finally:
                window.close()


if __name__=="__main__":
    unittest.main()
