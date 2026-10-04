from __future__ import annotations

from dataclasses import replace
from PySide6 import QtCore, QtGui, QtWidgets
import math
from ..api import load_preset
from ..config import RunConfig
from ..experiments import preset_names
from .instruments import InstrumentPanel
from .scene import Camera, scene_bounds
from .worker import SimulationWorker


class ApparatusView(QtWidgets.QWidget):
    def __init__(self, sim):
        super().__init__()
        shaft = getattr(sim.world.mechanism, "shaft", None)
        self.controlled_shaft = (shaft is not None and
                                 shaft.prescribed_omega is not None)
        self.bounds=scene_bounds(sim)
        state=sim.world.bodies
        self.shape=tuple(int(value) for value in state.shape)
        self.radius=tuple(float(value) for value in state.radius)
        self.polygons={i:tuple(tuple(float(v) for v in vertex)
                                for vertex in polygon)
                       for i,polygon in state.polygons.items()}
        self.walls=[]
        for wall in sim.world.walls:
            if hasattr(wall,"start"):
                start,end=tuple(wall.start),tuple(wall.end)
            else:
                nx,ny=wall.inward_normal
                tx,ty=-ny,nx
                x,y=wall.point_at(0.)
                half=wall.length/2
                start,end=(x-half*tx,y-half*ty),(x+half*tx,y+half*ty)
            self.walls.append((start,end,wall.kind_at(0.).value))
        self.snapshot=sim.snapshot()
        self.zoom=1.0
        self.pan_x=0.0
        self.pan_y=0.0
        self.selected_particle=None
        self._drag_at=None
        self.setMinimumSize(500,340)

    def set_snapshot(self, snapshot):
        self.snapshot=snapshot
        self.update()

    def _point(self, camera, point):
        return QtCore.QPointF(*camera.map(*point))

    def camera(self):
        return Camera(self.bounds,self.width(),self.height(),zoom=self.zoom,
                      pan_x=self.pan_x,pan_y=self.pan_y)

    def reset_camera(self):
        self.zoom=1.0
        self.pan_x=self.pan_y=0.0
        self.update()

    def wheelEvent(self,event):
        camera=self.camera()
        anchor=camera.unmap(event.position().x(),event.position().y())
        factor=1.2 if event.angleDelta().y()>0 else 1/1.2
        self.zoom=max(.5,min(12.,self.zoom*factor))
        moved=self.camera().unmap(event.position().x(),event.position().y())
        self.pan_x += anchor[0]-moved[0]
        self.pan_y += anchor[1]-moved[1]
        self.update()
        event.accept()

    def mousePressEvent(self,event):
        if event.button() != QtCore.Qt.MouseButton.LeftButton:
            return super().mousePressEvent(event)
        camera=self.camera()
        x,y=event.position().x(),event.position().y()
        distances=[math.hypot(x-camera.map(*point)[0],
                              y-camera.map(*point)[1])
                   for point in self.snapshot.position]
        nearest=min(range(len(distances)),key=distances.__getitem__) if distances else None
        if nearest is not None and distances[nearest] <= max(
                10.,self.radius[nearest]*camera.scale+5.):
            self.selected_particle=nearest
            self._drag_at=None
        else:
            self.selected_particle=None
            self._drag_at=event.position()
        self.update()
        event.accept()

    def mouseMoveEvent(self,event):
        if self._drag_at is None:
            return super().mouseMoveEvent(event)
        dx=event.position().x()-self._drag_at.x()
        dy=event.position().y()-self._drag_at.y()
        scale=self.camera().scale
        self.pan_x -= dx/scale
        self.pan_y += dy/scale
        self._drag_at=event.position()
        self.update()
        event.accept()

    def mouseReleaseEvent(self,event):
        self._drag_at=None
        super().mouseReleaseEvent(event)

    def component_visible(self, component):
        # The stored spring path is a phase marker, not a literal deforming
        # spring. Keep its energy in the instruments without drawing a coil.
        return component.name != "shaft_spring"

    def _component(self, painter, camera, component):
        colors = {
            "cylinder":"#8997a7", "piston":"#f0c46c",
            "piston_cam":"#75bdd5", "cam_follower":"#e8d8a8",
            "selector_cam":"#ab9edb", "selector_shoe":"#e8d8a8",
            "thermal_wall":{"hot":"#f2a65a", "cold":"#76bce3"}.get(
                component.state,"#d5dce4"),
            "flywheel":"#75bdd5" if self.controlled_shaft else "#d8dee7",
            "flywheel_spoke":"#f0c46c",
            "shaft_spring":"#9ed9a2", "load":"#dc9c9c",
        }
        color = QtGui.QColor(colors.get(component.name,"#d5dce4"))
        painter.setPen(QtGui.QPen(color,2.5 if component.name in
                                  ("piston","thermal_wall") else 1.6))
        painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
        points = [self._point(camera,point) for point in component.points]
        if component.kind == "point":
            painter.setBrush(color)
            painter.drawEllipse(points[0],5,5)
        elif len(points) >= 2:
            painter.drawPolyline(QtGui.QPolygonF(points))

    def _wall(self, painter, camera, wall):
        start,end,kind=wall
        color={"hot":"#f2a65a","cold":"#76bce3"}.get(kind,"#8997a7")
        painter.setPen(QtGui.QPen(QtGui.QColor(color),2))
        painter.drawLine(self._point(camera,start),self._point(camera,end))

    def paintEvent(self, event):
        painter=QtGui.QPainter(self)
        painter.fillRect(self.rect(),QtGui.QColor("#15191f"))
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        camera=self.camera()
        snap=self.snapshot
        if snap.apparatus:
            for component in snap.apparatus:
                if self.component_visible(component):
                    self._component(painter,camera,component)
        else:
            for wall in self.walls:
                self._wall(painter,camera,wall)
        painter.setBrush(QtGui.QColor("#e8b44d"))
        painter.setPen(QtGui.QPen(QtGui.QColor("#f3e7c3"),1))
        for i,(x,y) in enumerate(snap.position):
            center=self._point(camera,(x,y))
            if self.shape[i]==0:
                radius=self.radius[i]*camera.scale
                painter.drawEllipse(center,radius,radius)
            else:
                cosine,sine=math.cos(snap.angle[i]),math.sin(snap.angle[i])
                vertices=[]
                for vx,vy in self.polygons[i]:
                    vertices.append(self._point(camera,(x+cosine*vx-sine*vy,
                                                          y+sine*vx+cosine*vy)))
                painter.drawPolygon(QtGui.QPolygonF(vertices))
                painter.drawLine(center,vertices[0])
            if i == self.selected_particle:
                painter.setPen(QtGui.QPen(QtGui.QColor("#ffffff"),2))
                painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
                painter.drawEllipse(center,max(8.,self.radius[i]*camera.scale+6.),
                                    max(8.,self.radius[i]*camera.scale+6.))
                painter.setBrush(QtGui.QColor("#e8b44d"))
                painter.setPen(QtGui.QPen(QtGui.QColor("#f3e7c3"),1))
        if snap.apparatus:
            painter.setPen(QtGui.QColor("#d5dce4"))
            parts={part.name:part for part in snap.apparatus}
            def center(name):
                points=parts[name].points
                return ((min(x for x,_ in points)+max(x for x,_ in points))/2,
                        (min(y for _,y in points)+max(y for _,y in points))/2)
            shaft_center=parts["flywheel_spoke"].points[0]
            labels=[("Piston",parts["piston"].points[1],-50,-8),
                    ("Motor · fixed speed" if self.controlled_shaft else "Flywheel · free shaft",
                     shaft_center,-55,43),
                    ("Piston cam",center("piston_cam"),-36,43),
                    ("Thermal selector",center("selector_cam"),-51,43),
                    (f"Load W={parts['load'].work_output:.3g}",
                     parts["load"].points[-1],7,5)]
            if parts["shaft_spring"].state != "inactive":
                labels.append((f"Torsion spring E={parts['shaft_spring'].potential_energy:.3g}",
                               shaft_center,-65,-33))
            for label,point,offset_x,offset_y in labels:
                screen=self._point(camera,point)
                painter.drawText(screen+QtCore.QPointF(offset_x,offset_y),label)
            mode = "motor-driven" if self.controlled_shaft else "free shaft"
            painter.drawText(16,24,
                f"{snap.branch.replace('_',' ')} · {mode} · phase {snap.shaft_phase:.2f} rad")
        if self.selected_particle is not None:
            i=self.selected_particle
            painter.setPen(QtGui.QColor("#f3e7c3"))
            painter.drawText(16,self.height()-14,
                f"Particle {i} · v={math.hypot(*snap.velocity[i]):.3g} · "
                f"ω={snap.omega[i]:.3g}")
        painter.end()


class MainWindow(QtWidgets.QMainWindow):
    play_requested=QtCore.Signal(bool)
    step_requested=QtCore.Signal()
    advance_requested=QtCore.Signal(float)
    branch_requested=QtCore.Signal()
    playback_rate_requested=QtCore.Signal(float)
    shaft_speed_requested=QtCore.Signal(float)
    reset_requested=QtCore.Signal(int)
    configuration_requested=QtCore.Signal(object)
    frame_received=QtCore.Signal()
    stop_requested=QtCore.Signal()

    def __init__(self,preset,autoplay=True,transient_cycles=2,
                 efficiency_min_cycles=8,shaft_speed=.15,shaft_mode="controlled"):
        super().__init__(); self.setWindowTitle("Microscopic Thermodynamics Laboratory" +
            (" — experimental triangle CCD" if preset=="carnot_triangles" else ""))
        config=RunConfig(preset=preset,
                         particles=96 if preset=="carnot_triangles" else 48,
                         max_horizon=.02,shaft_speed=shaft_speed,
                         shaft_mode=shaft_mode)
        simulation=load_preset(config)
        self.config=config
        self._generation=0
        self.view=ApparatusView(simulation)
        self.instruments=InstrumentPanel()
        self.splitter=QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        self.splitter.addWidget(self.view)
        self.splitter.addWidget(self.instruments)
        self.splitter.setStretchFactor(0,2)
        self.splitter.setStretchFactor(1,1)
        self.splitter.setSizes([720,360])
        self.setCentralWidget(self.splitter)
        self.resize(1100,700)
        bar=self.addToolBar("Simulation")
        bar.addWidget(QtWidgets.QLabel(" Scene "))
        self.preset_choice=QtWidgets.QComboBox()
        for name in preset_names():
            self.preset_choice.addItem(name)
        self.preset_choice.setCurrentText(preset)
        self.preset_choice.currentTextChanged.connect(self.change_preset)
        bar.addWidget(self.preset_choice)
        self.play=bar.addAction("Pause" if autoplay else "Play")
        self.play.setCheckable(True)
        self.play.setChecked(not autoplay)
        self.play.toggled.connect(self._toggle_playback)
        step=bar.addAction("Step collision"); step.triggered.connect(self.step_collision)
        self.duration_step=QtWidgets.QDoubleSpinBox()
        self.duration_step.setRange(.001,100.)
        self.duration_step.setDecimals(3)
        self.duration_step.setSingleStep(.01)
        self.duration_step.setValue(.1)
        self.duration_step.setSuffix(" s")
        bar.addWidget(self.duration_step)
        fixed=bar.addAction("Step duration")
        fixed.triggered.connect(self.step_duration)
        self.branch_action=bar.addAction("Step branch")
        self.branch_action.setEnabled(preset.startswith("carnot_") and
                                      shaft_mode=="controlled")
        self.branch_action.triggered.connect(self.step_branch)
        bar.addWidget(QtWidgets.QLabel(" Playback "))
        self.playback_rate=QtWidgets.QComboBox()
        for rate in (.25,.5,1.,2.,4.):
            self.playback_rate.addItem(f"{rate:g}×",rate)
        self.playback_rate.setCurrentIndex(2)
        self.playback_rate.currentIndexChanged.connect(
            lambda _: self.playback_rate_requested.emit(
                float(self.playback_rate.currentData())))
        bar.addWidget(self.playback_rate)
        bar.addWidget(QtWidgets.QLabel(" Shaft "))
        self.shaft_speed=QtWidgets.QDoubleSpinBox()
        self.shaft_speed.setRange(.005,3.)
        self.shaft_speed.setDecimals(3)
        self.shaft_speed.setSingleStep(.025)
        self.shaft_speed.setValue(shaft_speed)
        self.shaft_speed.setSuffix(" rad/s")
        self.shaft_speed.setEnabled(preset.startswith("carnot_") and
                                    shaft_mode=="controlled")
        self.shaft_speed.setToolTip(
            "Changes the physical experiment; recorded as a motor intervention.")
        self.shaft_speed.editingFinished.connect(self.change_shaft_speed)
        bar.addWidget(self.shaft_speed)
        bar.addWidget(QtWidgets.QLabel(" Seed "))
        self.seed=QtWidgets.QSpinBox()
        self.seed.setRange(0,999999999)
        self.seed.setValue(config.seed)
        bar.addWidget(self.seed)
        reset=bar.addAction("Reset to seed")
        reset.triggered.connect(self.reset_to_seed)
        camera_reset=bar.addAction("Fit scene")
        camera_reset.triggered.connect(lambda: self.view.reset_camera())
        self.status=self.statusBar()
        self.thread=QtCore.QThread(self)
        self.worker=SimulationWorker(
            simulation,autoplay=autoplay,transient_cycles=transient_cycles,
            efficiency_min_cycles=efficiency_min_cycles,config=config)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.start)
        self.play_requested.connect(self.worker.set_playing)
        self.step_requested.connect(self.worker.step_collision)
        self.advance_requested.connect(self.worker.advance_duration)
        self.branch_requested.connect(self.worker.step_branch)
        self.playback_rate_requested.connect(self.worker.set_playback_rate)
        self.shaft_speed_requested.connect(self.worker.set_shaft_speed)
        self.reset_requested.connect(self.worker.reset_seed)
        self.configuration_requested.connect(self.worker.reset_configuration)
        self.frame_received.connect(self.worker.frame_received)
        self.stop_requested.connect(self.worker.stop)
        self.worker.snapshot_ready.connect(self.on_snapshot)
        self.worker.failed.connect(self.on_failure)
        self.worker.finished.connect(self.thread.quit,
                                     QtCore.Qt.ConnectionType.DirectConnection)
        self.thread.finished.connect(self.worker.deleteLater)
        self.thread.start()

    @QtCore.Slot(object)
    def on_snapshot(self, frame):
        if frame.generation != self._generation:
            self.frame_received.emit()
            return
        snapshot=frame.snapshot
        self.view.set_snapshot(snapshot)
        self.instruments.set_frame(frame.instruments)
        self.instruments.set_rates(frame.rates)
        self.status.showMessage(
            f"{'PAUSED' if self.play.isChecked() else 'RUNNING'}  "
            f"t={snapshot.time:.3f}  events={snapshot.event_count}  "
            f"T={snapshot.translational_temperature:.3f}  "
            f"R_E={snapshot.energy_residual:.2e}")
        self.frame_received.emit()

    @QtCore.Slot(str,str)
    def on_failure(self, message, path):
        self.play.setChecked(True)
        detail=f"; diagnostic: {path}" if path else ""
        self.status.showMessage(f"PAUSED: {message}{detail}")
        QtWidgets.QMessageBox.critical(self,"Numerical failure",
            f"The scientific run was paused.\n{message}\n"
            f"Diagnostic and checkpoint: {path}" if path else message)

    def step_collision(self):
        self.play.setChecked(True)
        self.step_requested.emit()

    def _toggle_playback(self, paused):
        self.play.setText("Play" if paused else "Pause")
        self.play_requested.emit(not paused)

    def step_duration(self):
        self.play.setChecked(True)
        self.advance_requested.emit(self.duration_step.value())

    def step_branch(self):
        self.play.setChecked(True)
        self.branch_requested.emit()

    def change_shaft_speed(self):
        if not self.shaft_speed.isEnabled():
            return
        self.play.setChecked(True)
        self.config=replace(self.config,shaft_speed=self.shaft_speed.value())
        self.shaft_speed_requested.emit(self.shaft_speed.value())

    def reset_to_seed(self):
        self.play.setChecked(True)
        self.view.selected_particle=None
        self.config=replace(self.config,seed=self.seed.value())
        self._generation += 1
        self.reset_requested.emit(self.seed.value())

    def change_preset(self, preset):
        if preset == self.config.preset:
            return
        self.play.setChecked(True)
        config=RunConfig(preset=preset,seed=self.seed.value(),
                         particles=96 if preset=="carnot_triangles" else 48,
                         max_horizon=.02,shaft_speed=self.shaft_speed.value(),
                         shaft_mode=self.config.shaft_mode)
        replacement=ApparatusView(load_preset(config))
        previous=self.splitter.replaceWidget(0,replacement)
        self.view=replacement
        previous.deleteLater()
        self.config=config
        can_drive=preset.startswith("carnot_") and config.shaft_mode=="controlled"
        self.branch_action.setEnabled(can_drive)
        self.shaft_speed.setEnabled(can_drive)
        self.setWindowTitle("Microscopic Thermodynamics Laboratory"+
            (" — experimental triangle CCD" if preset=="carnot_triangles" else ""))
        self._generation += 1
        self.configuration_requested.emit(config)

    def closeEvent(self, event):
        self.stop_requested.emit()
        self.thread.wait()
        super().closeEvent(event)


def launch(preset="gas_box",transient_cycles=2,efficiency_min_cycles=8,
           shaft_speed=.15,shaft_mode="controlled"):
    app=QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window=MainWindow(preset,transient_cycles=transient_cycles,
                      efficiency_min_cycles=efficiency_min_cycles,
                      shaft_speed=shaft_speed,shaft_mode=shaft_mode)
    window.show(); return app.exec()
