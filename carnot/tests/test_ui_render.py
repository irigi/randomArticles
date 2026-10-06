"""Offscreen smoke check for the optional Qt apparatus renderer."""

import os
import time
import math
import unittest
import tempfile
from pathlib import Path
import numpy as np

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
try:
    from PySide6 import QtCore, QtTest, QtWidgets
except ImportError:
    QtWidgets = None

from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.ui.scene import Camera, scene_bounds


@unittest.skipIf(QtWidgets is None, "PySide6 is not installed")
class UiRenderTests(unittest.TestCase):
    def wait_until(self, app, condition, timeout=5.):
        deadline=time.monotonic()+timeout
        while time.monotonic()<deadline:
            app.processEvents()
            if condition():
                return
            time.sleep(.005)
        self.fail("timed out waiting for Qt worker snapshot")

    def test_small_carnot_view_draws_only_the_cylinder(self):
        from microthermo.ui.main_window import ApparatusView

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        sim = load_preset(RunConfig(preset="carnot_triangles",particles=8,seed=123,
                                    hot_jacket=True,cold_jacket=True))
        view = ApparatusView(sim)
        view.resize(500,200)
        view.show()
        app.processEvents()
        image = view.grab().toImage()
        self.assertEqual((image.width(),image.height()),(500,200))
        camera = Camera(scene_bounds(sim,cylinder_only=True),500,200)
        parts = {part.name:part for part in sim.snapshot().apparatus}
        self.assertTrue(view.controlled_shaft)
        self.assertTrue(view.cylinder_only)
        for name in ("piston_cam","selector_cam","flywheel","load","shaft_spring"):
            self.assertFalse(view.component_visible(parts[name]),name)
        # Piston, thermal wall, and the active hot jacket on the bottom wall.
        for name,point in (("piston",parts["piston"].points[0]),
                           ("thermal_wall",parts["thermal_wall"].points[0]),
                           ("jacket",(parts["piston"].points[0][0]/2,0.))):
            x,y = camera.map(*point)
            pixels = [image.pixelColor(int(x)+dx,int(y)+dy)
                      for dx in (-2,-1,0,1,2) for dy in (-2,-1,0,1,2)]
            self.assertTrue(any(color.red()+color.green()+color.blue() > 300
                                for color in pixels),name)
        view.close()

    def test_cycle_dial_tracks_sector_progress(self):
        from microthermo.ui.lab_view import CycleDial

        QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        sim = load_preset(RunConfig(preset="carnot_discs",particles=4,
                                    cam_fractions=(.4,.2,.2,.2)))
        dial = CycleDial(sim)
        dial.set_snapshot(sim.advance_to(.2*2*math.pi/.15))
        self.assertEqual(dial.branch,"hot")
        self.assertAlmostEqual(dial.phase_fraction(),.5,places=9)
        self.assertEqual([s[2] for s in dial.sectors],
                         ["hot","adiabatic_expansion","cold","adiabatic_compression"])
        self.assertGreater(dial.grab().toImage().width(),0)

    def test_camera_zoom_pan_and_inverse_mapping(self):
        from microthermo.ui.scene import Bounds

        camera=Camera(Bounds(-1.,-2.,3.,2.),800,600,zoom=2.,
                      pan_x=.25,pan_y=-.5)
        for point in ((0.,0.),(2.,1.),(-.5,-1.5)):
            np.testing.assert_allclose(camera.unmap(*camera.map(*point)),
                                       point,atol=1e-12)
        self.assertGreater(camera.scale,
                           Camera(camera.bounds,800,600).scale)

    def test_apparatus_and_particle_inspection(self):
        from microthermo.ui.main_window import MainWindow

        app=QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        window=MainWindow("carnot_discs",autoplay=False)
        try:
            window.show()
            self.wait_until(app,lambda: window.readout.values["Area"].text()!="—")
            self.assertFalse(window.readout.inspector.isVisibleTo(window))
            piston=next(part for part in window.view.snapshot.apparatus
                        if part.name=="piston")
            x,y=window.view.camera().map(*piston.points[0])
            QtTest.QTest.mouseClick(window.view,QtCore.Qt.MouseButton.LeftButton,
                                    pos=QtCore.QPoint(round(x),round(y)-20))
            self.assertIn("piston",window.readout.inspector.text().lower())
            self.assertTrue(window.readout.inspector.isVisibleTo(window))
            x,y=window.view.camera().map(*window.view.snapshot.position[0])
            QtTest.QTest.mouseClick(window.view,QtCore.Qt.MouseButton.LeftButton,
                                    pos=QtCore.QPoint(round(x),round(y)))
            self.assertIn("Particle 0",window.readout.inspector.text())
            self.assertIn("velocity",window.readout.inspector.text())
        finally:
            window.close()

    def test_dropped_render_frames_preserve_physics_horizon(self):
        from microthermo.ui.worker import SimulationWorker

        QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        config=RunConfig(preset="gas_box",particles=8,seed=124,max_horizon=.02)
        simulation=load_preset(config)
        worker=SimulationWorker(simulation,autoplay=False,config=config)
        frames=[]
        worker.snapshot_ready.connect(frames.append)
        worker.start()
        worker.advance_duration(.01)
        worker.advance_duration(.02)
        self.assertEqual(len(frames),1)
        worker.frame_received()
        self.assertEqual(len(frames),2)
        reference=load_preset(config)
        expected=reference.advance_to(.03)
        actual=frames[-1].snapshot
        np.testing.assert_array_equal(actual.position,expected.position)
        np.testing.assert_array_equal(actual.velocity,expected.velocity)
        self.assertEqual(actual.event_count,expected.event_count)
        self.assertAlmostEqual(actual.time,.03,places=12)
        worker.stop()

    def test_free_shaft_and_active_spring_are_labeled_by_state(self):
        from microthermo.ui.main_window import ApparatusView

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        sim = load_preset(RunConfig(preset="carnot_discs", particles=8,
                                    shaft_mode="free"))
        sim.world.mechanism.shaft.spring_k = 1.0
        view = ApparatusView(sim)
        parts = {part.name: part for part in sim.snapshot().apparatus}
        self.assertFalse(view.controlled_shaft)
        self.assertEqual(parts["shaft_spring"].state, "torsion")
        self.assertFalse(view.component_visible(parts["shaft_spring"]))
        view.close()

    def test_worker_command_matches_headless_and_snapshot_is_immutable(self):
        from microthermo.ui.main_window import MainWindow

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        window=MainWindow("carnot_triangles",autoplay=False)
        try:
            window.show()
            self.wait_until(app,lambda: window.thread.isRunning())
            window.advance_requested.emit(.05)
            self.wait_until(app,lambda: window.view.snapshot.time>=.05)
            snap=window.view.snapshot
            reference=load_preset(RunConfig(preset="carnot_triangles",
                                                particles=96,max_horizon=.02))
            expected=reference.advance_to(.05)
            np.testing.assert_array_equal(snap.position,expected.position)
            np.testing.assert_array_equal(snap.velocity,expected.velocity)
            self.assertEqual(snap.event_count,expected.event_count)
            self.assertEqual(snap.shaft_phase,expected.shaft_phase)
            with self.assertRaises(ValueError):
                snap.position[0,0]=0.
            with self.assertRaises(ValueError):
                snap.position.setflags(write=True)
            self.assertIs(window.worker.thread(),window.thread)
        finally:
            window.close()
        self.assertFalse(window.thread.isRunning())

    def test_worker_pause_resume_and_explicit_step(self):
        from microthermo.ui.main_window import MainWindow

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        window=MainWindow("gas_box",autoplay=False)
        try:
            window.show()
            self.wait_until(app,lambda: window.thread.isRunning())
            window.play_requested.emit(True)
            self.wait_until(app,lambda: window.view.snapshot.time>=.02)
            self.wait_until(app,lambda:
                window.readout.diagnostic_values["Achieved playback"].text()
                .endswith("×"))
            self.assertIn("events/s",window.readout.diagnostic_values[
                "Event throughput"].text())
            window.play_requested.emit(False)
            self.wait_until(app,lambda:
                window.readout.diagnostic_values["Achieved playback"].text()
                == "Paused")
            paused=window.view.snapshot.time
            time.sleep(.05)
            app.processEvents()
            self.assertEqual(window.view.snapshot.time,paused)
            self.assertEqual(window.readout.diagnostic_values[
                "Achieved playback"].text(),"Paused")
            window.advance_requested.emit(.01)
            self.wait_until(app,lambda: window.view.snapshot.time>paused)
            self.assertAlmostEqual(window.view.snapshot.time,paused+.01,places=12)
        finally:
            window.close()

    def test_worker_drives_the_edmd_engine(self):
        from microthermo.core.edmd import numba_available
        from microthermo.ui.main_window import MainWindow

        if not numba_available():
            self.skipTest("Numba optional dependency not installed")
        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        window=MainWindow("gas_box",autoplay=False,
                          config=RunConfig(preset="gas_box",particles=48,engine="edmd"))
        try:
            window.show()
            self.wait_until(app,lambda: window.thread.isRunning())
            window.play_requested.emit(True)
            self.wait_until(app,lambda: window.view.snapshot.time>=.5)
            window.play_requested.emit(False)
            self.wait_until(app,lambda:
                window.readout.diagnostic_values["Achieved playback"].text()
                == "Paused")
            paused=window.view.snapshot.time
            window.advance_requested.emit(.01)
            self.wait_until(app,lambda: window.view.snapshot.time>paused)
            self.assertAlmostEqual(window.view.snapshot.time,paused+.01,places=12)
            self.assertGreater(window.view.snapshot.event_count,0)
        finally:
            window.close()

    def test_desktop_duration_rate_reset_and_physical_speed_are_distinct(self):
        from microthermo.ui.main_window import MainWindow

        app=QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        window=MainWindow("carnot_discs",autoplay=False)
        try:
            window.show()
            self.wait_until(app,lambda: window.view.snapshot.time==0 and
                            window.readout.values["Area"].text() != "—")
            window.playback_rate.setCurrentIndex(window.playback_rate.findData(4.))
            self.wait_until(app,lambda: abs(window.worker._physical_step-.064)<1e-12)
            self.assertAlmostEqual(window.view.snapshot.shaft_phase,0.)
            window.duration_step.setValue(.05)
            window.step_duration()
            self.wait_until(app,lambda: window.view.snapshot.time>=.05)
            self.assertAlmostEqual(window.view.snapshot.time,.05,places=12)
            window.shaft_speed.setValue(.3)
            window.change_shaft_speed()
            self.wait_until(app,lambda: window.view.snapshot.event_count>=1)
            self.assertEqual(window.worker._simulation.events[-1].kind,
                             "shaft_speed_change")
            self.assertLess(abs(window.view.snapshot.energy_residual),1e-10)
            window.seed.setValue(124)
            window.reset_to_seed()
            self.wait_until(app,lambda: window._generation==1 and
                            window.view.snapshot.time==0 and
                            window.worker._generation==1)
            expected=load_preset(RunConfig(preset="carnot_discs",particles=48,
                                                seed=124,max_horizon=.02,
                                                shaft_speed=.3)).snapshot()
            np.testing.assert_array_equal(window.view.snapshot.position,
                                          expected.position)
            self.assertEqual(window.status.currentMessage().split()[0],"PAUSED")
            window.preset_choice.setCurrentText("gas_box")
            self.wait_until(app,lambda: window._generation==2 and
                            window.worker._generation==2 and
                            window.view.snapshot.time==0)
            gas=load_preset(RunConfig(preset="gas_box",particles=48,
                                           seed=124,max_horizon=.02)).snapshot()
            np.testing.assert_array_equal(window.view.snapshot.position,
                                          gas.position)
            self.assertFalse(window.shaft_speed.isEnabled())
            self.assertFalse(window.branch_action.isEnabled())
        finally:
            window.close()

    def test_branch_step_stops_at_next_controlled_boundary(self):
        from microthermo.ui.worker import SimulationWorker

        QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        config=RunConfig(preset="carnot_discs",particles=2,seed=123,
                         shaft_speed=.15,max_horizon=.05)
        sim=load_preset(config)
        worker=SimulationWorker(sim,autoplay=False,config=config)
        worker.step_branch()
        self.assertAlmostEqual(sim.time,math.pi/(2*.15),delta=1e-9)
        self.assertEqual(sim.snapshot().branch,"adiabatic_expansion")

    def test_desktop_checkpoint_config_and_exports_round_trip(self):
        from microthermo.ui.main_window import MainWindow

        app=QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        window=MainWindow("carnot_discs",autoplay=False)
        try:
            window.show()
            self.wait_until(app,lambda: window.readout.values["Area"].text() != "—")
            window.advance_requested.emit(.05)
            self.wait_until(app,lambda: window.view.snapshot.time>=.05)
            original=window.view.snapshot
            with tempfile.TemporaryDirectory() as folder:
                root=Path(folder)
                window.save_configuration(root/"config.json")
                window.save_checkpoint(root/"checkpoint.pkl.gz")
                self.wait_until(app,lambda: (root/"checkpoint.pkl.gz").exists())
                window.export_run(root/"export")
                self.wait_until(app,lambda: (root/"export"/"summary.json").exists())
                window.export_plots(root/"plots.png")
                self.assertTrue((root/"plots_pressure_area.png").exists())
                self.assertTrue((root/"plots_temperatures.png").exists())
                self.assertTrue((root/"plots_temperature_entropy.png").exists())
                window.start_recording(root/"recording")
                window._capture_recording_frame()
                window.stop_recording()
                self.assertEqual(window.state_label.text(),"PAUSED")
                import json
                manifest=json.loads((root/"recording"/"recording.json").read_text())
                self.assertEqual(len(manifest["frames"]),2)
                self.assertEqual(manifest["frames"][0]["physical_time"],.05)
                window.advance_requested.emit(.02)
                self.wait_until(app,lambda: window.view.snapshot.time>=.07)
                window.open_checkpoint(root/"checkpoint.pkl.gz")
                self.wait_until(app,lambda: window._generation==1 and
                                window.view.snapshot.time==original.time)
                np.testing.assert_array_equal(window.view.snapshot.position,
                                              original.position)
                self.assertEqual(window.view.snapshot.event_count,
                                 original.event_count)
                window.open_configuration(root/"config.json")
                self.wait_until(app,lambda: window._generation==2 and
                                window.view.snapshot.time==0)
        finally:
            window.close()

    def test_lab_layout_fits_a_laptop_screen(self):
        from microthermo.ui.main_window import MainWindow

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        window=MainWindow("carnot_triangles",autoplay=False)
        try:
            window.resize(1280,800)
            window.show()
            self.wait_until(app,lambda: window.readout.values["Area"].text() != "—")
            self.assertGreaterEqual(window.view.width(),300)
            self.assertGreaterEqual(window.readout.width(),320)
            for toolbar in window.findChildren(QtWidgets.QToolBar):
                self.assertTrue(all(toolbar.actionGeometry(action).right()
                                    <= toolbar.width()
                                    for action in toolbar.actions()))
            self.assertTrue(window.view.isVisible())
            self.assertTrue(window.lab.dial.isVisible())
            self.assertTrue(window.diagrams.isVisible())
            self.assertTrue(window.readout.isVisible())
            window.advance_requested.emit(.2)
            self.wait_until(app,lambda: window.view.snapshot.time>=.2)
            def temperatures_drawn():
                x=window.readout.temp_trans.getData()[0]
                return x is not None and len(x)>=2
            self.wait_until(app,temperatures_drawn)
            area_x,pressure_y=window.diagrams.pv_lines["hot"].getData()
            entropy_x,temperature_y=window.diagrams.ts_lines["hot"].getData()
            self.assertGreaterEqual(len(area_x),2)
            self.assertEqual(len(area_x),len(pressure_y))
            self.assertEqual(len(entropy_x),len(temperature_y))
            self.assertEqual(len(window.diagrams.pv_reference["hot"].getData()[0]),65)
            # The ideal T–S rectangle starts at zero entropy on the hot isotherm.
            ts_x,ts_y=window.diagrams.ts_reference["hot"].getData()
            self.assertEqual(ts_x[0],0.)
            self.assertEqual(ts_y[0],1.5)
            self.assertIn("Q_H",window.diagrams.current_cycle.text())
            self.assertNotEqual(window.readout.values["Gas energy"].text(),"—")
            self.assertNotEqual(window.readout.values["First-law residual"].text(),"—")
            self.assertFalse(window.readout.diagnostics.body.isVisible())
            self.assertFalse(window.readout.mechanism.body.isVisible())
            window.readout.diagnostics_toggle.setChecked(True)
            self.assertTrue(window.readout.diagnostics.body.isVisible())
            self.assertNotEqual(window.readout.diagnostic_values[
                "Max penetration"].text(),"—")
            # The right column fits without scrolling at this size.
            scroll=window.readout.findChild(QtWidgets.QScrollArea)
            window.readout.diagnostics_toggle.setChecked(False)
            app.processEvents()
            self.assertEqual(scroll.verticalScrollBar().maximum(),0)
            image=window.grab().toImage()
            self.assertEqual((image.width(),image.height()),
                             (window.width(),window.height()))
        finally:
            window.close()

    def test_non_carnot_scene_hides_cycle_diagrams(self):
        from microthermo.ui.main_window import MainWindow

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        window=MainWindow("gas_box",autoplay=False)
        try:
            window.show()
            self.wait_until(app,lambda: window.readout.values["Gas energy"].text() != "—")
            self.assertFalse(window.diagrams.isVisible())
            self.assertIsNone(window.lab.dial)
        finally:
            window.close()

    def test_completed_cycle_totals_reach_readout_and_scoreboard(self):
        from microthermo.measurements.live import LiveInstruments
        from microthermo.ui.instruments import (CycleDiagrams, InstrumentSeries,
                                                ReadoutPanel)

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        sim=load_preset(RunConfig(preset="carnot_discs",particles=2,
                                      seed=123,max_horizon=.1))
        snapshot=sim.advance_to(2*math.pi/.15)
        frame=LiveInstruments().observe(sim,snapshot)
        panel=ReadoutPanel()
        panel.set_sample(frame.current)
        self.assertTrue(panel.values["Time / cycles"].text().endswith("/ 1"))
        self.assertNotEqual(panel.values["Q_H / Q_C"].text(),"—")
        self.assertNotEqual(panel.diagnostic_values["Last cycle gas residual"].text(),"—")
        self.assertIn("±",panel.values["Work by gas: ledger / ∮P dA"].text())
        panel.close()
        diagrams=CycleDiagrams()
        diagrams.set_reference(frame.ideal_reference)
        diagrams.update_series(InstrumentSeries.from_history(frame.history),
                               frame.current,frame.particles,frame.degrees_of_freedom)
        diagrams.update_scoreboard(frame.current,frame.cycles,frame.cycle_start,2,.5)
        self.assertEqual(diagrams.table.rowCount(),1)
        # Engine sign: the table shows work delivered by the gas.
        cycle=frame.cycles[0]
        self.assertEqual(diagrams.table.item(0,3).text(),
                         f"{-cycle.piston_work_on_gas:.4g}")
        self.assertIn("warm-up",diagrams.table.item(0,0).text())
        diagrams.close()


if __name__ == "__main__":
    unittest.main()
