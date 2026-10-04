"""Offscreen smoke check for the optional Qt apparatus renderer."""

import os
import time
import math
import unittest
import numpy as np

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
try:
    from PySide6 import QtWidgets
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

    def test_small_carnot_view_draws_cylinder_and_drive(self):
        from microthermo.ui.main_window import ApparatusView

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        sim = load_preset(RunConfig(preset="carnot_triangles",particles=8,seed=123))
        view = ApparatusView(sim)
        view.resize(500,340)
        view.show()
        app.processEvents()
        image = view.grab().toImage()
        self.assertEqual((image.width(),image.height()),(500,340))
        camera = Camera(scene_bounds(sim),500,340)
        parts = {part.name:part for part in sim.snapshot().apparatus}
        self.assertTrue(view.controlled_shaft)
        self.assertFalse(view.component_visible(parts["shaft_spring"]))
        for name,index in (("piston",0),("thermal_wall",0),
                           ("piston_cam",20),("selector_cam",20),
                           ("flywheel",20),
                           ("load",1)):
            x,y = camera.map(*parts[name].points[index])
            pixels = [image.pixelColor(int(x)+dx,int(y)+dy)
                      for dx in (-2,-1,0,1,2) for dy in (-2,-1,0,1,2)]
            self.assertTrue(any(color.red()+color.green()+color.blue() > 300
                                for color in pixels),name)
        view.close()

    def test_camera_zoom_pan_and_inverse_mapping(self):
        from microthermo.ui.scene import Bounds

        camera=Camera(Bounds(-1.,-2.,3.,2.),800,600,zoom=2.,
                      pan_x=.25,pan_y=-.5)
        for point in ((0.,0.),(2.,1.),(-.5,-1.5)):
            np.testing.assert_allclose(camera.unmap(*camera.map(*point)),
                                       point,atol=1e-12)
        self.assertGreater(camera.scale,
                           Camera(camera.bounds,800,600).scale)

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
                window.instruments.diagnostic_values["Achieved playback"].text()
                .endswith("×"))
            self.assertIn("events/s",window.instruments.diagnostic_values[
                "Event throughput"].text())
            window.play_requested.emit(False)
            self.wait_until(app,lambda:
                window.instruments.diagnostic_values["Achieved playback"].text()
                == "Paused")
            paused=window.view.snapshot.time
            time.sleep(.05)
            app.processEvents()
            self.assertEqual(window.view.snapshot.time,paused)
            self.assertEqual(window.instruments.diagnostic_values[
                "Achieved playback"].text(),"Paused")
            window.advance_requested.emit(.01)
            self.wait_until(app,lambda: window.view.snapshot.time>paused)
            self.assertAlmostEqual(window.view.snapshot.time,paused+.01,places=12)
        finally:
            window.close()

    def test_desktop_duration_rate_reset_and_physical_speed_are_distinct(self):
        from microthermo.ui.main_window import MainWindow

        app=QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        window=MainWindow("carnot_discs",autoplay=False)
        try:
            window.show()
            self.wait_until(app,lambda: window.view.snapshot.time==0 and
                            window.instruments.values["Area"].text() != "—")
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

    def test_splitter_and_live_plots_fit_minimum_window(self):
        from microthermo.ui.main_window import MainWindow

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        window=MainWindow("carnot_triangles",autoplay=False)
        try:
            window.resize(850,520)
            window.show()
            self.wait_until(app,lambda: window.instruments.values["Area"].text() != "—")
            self.assertGreaterEqual(window.view.width(),500)
            self.assertGreaterEqual(window.instruments.width(),320)
            self.assertTrue(window.view.isVisible())
            self.assertTrue(window.instruments.isVisible())
            window.advance_requested.emit(.05)
            self.wait_until(app,lambda: window.view.snapshot.time>=.05)
            area_x,pressure_y=window.instruments.pa_line.getData()
            temp_x,temp_y=window.instruments.temp_trans.getData()
            self.assertGreaterEqual(len(area_x),1)
            self.assertGreaterEqual(len(temp_x),2)
            self.assertEqual(len(area_x),len(pressure_y))
            self.assertEqual(len(temp_x),len(temp_y))
            self.assertEqual(len(window.instruments.pa_reference["hot"].getData()[0]),65)
            self.assertIn("analytical point gas",window.instruments.pa_legend.text())
            self.assertNotEqual(window.instruments.values["Gas energy"].text(),"—")
            self.assertNotEqual(window.instruments.values["First-law residual"].text(),"—")
            self.assertFalse(window.instruments.diagnostics.isVisible())
            window.instruments.diagnostics_toggle.setChecked(True)
            self.assertTrue(window.instruments.diagnostics.isVisible())
            self.assertNotEqual(window.instruments.diagnostic_values[
                "Max penetration"].text(),"—")
            self.assertNotEqual(window.instruments.diagnostic_values[
                "Max cluster residual"].text(),"—")
            image=window.grab().toImage()
            self.assertEqual((image.width(),image.height()),
                             (window.width(),window.height()))
        finally:
            window.close()

    def test_completed_cycle_totals_reach_instrument_panel(self):
        from microthermo.measurements.live import LiveInstruments
        from microthermo.ui.instruments import InstrumentPanel

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        sim=load_preset(RunConfig(preset="carnot_discs",particles=2,
                                      seed=123,max_horizon=.1))
        snapshot=sim.advance_to(2*math.pi/.15)
        panel=InstrumentPanel()
        panel.set_frame(LiveInstruments().observe(sim,snapshot))
        self.assertEqual(panel.values["Completed cycles"].text(),"1")
        self.assertNotEqual(panel.values["Last cycle QH / QC"].text(),"—")
        self.assertNotEqual(panel.values["Last cycle gas residual"].text(),"—")
        self.assertNotEqual(panel.values["Cycle ∮P dA / smoothing bound"].text(),"—")
        panel.close()


if __name__ == "__main__":
    unittest.main()
