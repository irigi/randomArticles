from __future__ import annotations

from dataclasses import asdict, fields, replace
from concurrent.futures import ThreadPoolExecutor
from PySide6 import QtCore, QtGui, QtWidgets
import json
import math
from pathlib import Path
import pyqtgraph.exporters
import time
from ..api import load_preset
from ..config import RunConfig
from ..experiments import preset_names
from .instruments import InstrumentSeries
from .scene import Camera, scene_bounds
from .worker import SimulationWorker


class ApparatusView(QtWidgets.QWidget):
    inspected=QtCore.Signal(str)

    def __init__(self, sim):
        super().__init__()
        shaft = getattr(sim.world.mechanism, "shaft", None)
        self.controlled_shaft = (shaft is not None and
                                 shaft.prescribed_omega is not None)
        mechanism = sim.world.mechanism
        cam = getattr(mechanism, "cam", None)
        # Carnot scenes draw only the cylinder; the drive train is summarized
        # by the cycle dial. Its energies remain in the readout panel.
        self.cylinder_only = cam is not None and hasattr(mechanism, "apparatus_state")
        self.bounds=scene_bounds(sim, cylinder_only=self.cylinder_only)
        if self.cylinder_only:
            self.barrel = (max(cam.areas)/cam.height, cam.height)
            metadata = sim.world.metadata
            self.jackets = {"hot": bool(metadata.get("hot_jacket")),
                            "cold": bool(metadata.get("cold_jacket"))}
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
        self.selected_component=None
        self._drag_at=None
        self.setMinimumSize(*((300,120) if self.cylinder_only else (500,340)))

    def set_snapshot(self, snapshot):
        self.snapshot=snapshot
        if self.selected_particle is not None or self.selected_component is not None:
            self.inspected.emit(self._inspection_text())
        self.update()

    def _inspection_text(self):
        if self.selected_particle is not None:
            i=self.selected_particle
            x,y=self.snapshot.position[i]
            vx,vy=self.snapshot.velocity[i]
            return (f"Particle {i}: position ({x:.3g}, {y:.3g}), "
                    f"velocity ({vx:.3g}, {vy:.3g}), "
                    f"rotation {self.snapshot.omega[i]:.3g} rad/s. "
                    "Contacts exchange momentum and may exchange heat at "
                    "thermal walls.")
        part=next((p for p in self.snapshot.apparatus
                   if p.name==self.selected_component),None)
        if part is None:
            return "Click a particle or apparatus part to inspect it."
        explanation={
            "cylinder":"The gas volume changes as the piston moves.",
            "piston":"The piston transfers mechanical work to or from the gas.",
            "piston_cam":"The cam prescribes piston position from shaft phase.",
            "cam_follower":"The follower transmits the cam's prescribed motion.",
            "selector_cam":"The selector cam chooses the thermal contact sector.",
            "selector_shoe":"The shoe follows the active thermal sector.",
            "thermal_wall":"The active wall exchanges heat with the selected bath.",
            "flywheel":"The shaft stores rotational energy when running free.",
            "flywheel_spoke":"The spoke shows the current shaft phase.",
            "load":"The load receives mechanical output from shaft torque.",
        }.get(part.name,"")
        values=[]
        if part.state:
            values.append(f"state {part.state.replace('_',' ')}")
        if part.kinetic_energy:
            values.append(f"kinetic energy {part.kinetic_energy:.3g}")
        if part.potential_energy:
            values.append(f"potential energy {part.potential_energy:.3g}")
        if part.work_output:
            values.append(f"load output {part.work_output:.3g}")
        if part.angular_velocity:
            values.append(f"angular speed {part.angular_velocity:.3g} rad/s")
        return (part.name.replace("_"," ").title()+": "+explanation+
                (" Current "+", ".join(values)+"." if values else ""))

    def _nearest_component(self,x,y,camera):
        best=None
        for part in self.snapshot.apparatus:
            if not self.component_visible(part):
                continue
            points=[camera.map(*p) for p in part.points]
            distances=[]
            if len(points)==1:
                distances.append(math.hypot(x-points[0][0],y-points[0][1]))
            for (ax,ay),(bx,by) in zip(points[:-1],points[1:]):
                dx,dy=bx-ax,by-ay
                length=dx*dx+dy*dy
                fraction=max(0.,min(1.,((x-ax)*dx+(y-ay)*dy)/length)) if length else 0.
                distances.append(math.hypot(x-ax-fraction*dx,
                                            y-ay-fraction*dy))
            if distances:
                distance=min(distances)
                # The cylinder outline shares the piston and thermal wall.
                rank=(distance,part.name=="cylinder")
                if distance<=12. and (best is None or rank<best[0]):
                    best=(rank,part.name)
        return best[1] if best else None

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
            self.selected_component=None
            self._drag_at=None
        else:
            self.selected_particle=None
            self.selected_component=self._nearest_component(x,y,camera)
            self._drag_at=(event.position() if self.selected_component is None
                           else None)
        self.inspected.emit(self._inspection_text())
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
        if self.cylinder_only:
            return component.name in ("cylinder", "piston", "thermal_wall")
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

    def _paint_jackets(self, painter, camera):
        # Thermal jackets make the top and bottom walls exchange heat during
        # their sector; colour them while they are active.
        branch=self.snapshot.branch
        if branch not in ("hot","cold") or not self.jackets.get(branch):
            return
        piston=next(p for p in self.snapshot.apparatus if p.name=="piston")
        right=piston.points[0][0]
        top=max(y for _,y in piston.points)
        color=QtGui.QColor("#f2a65a" if branch=="hot" else "#76bce3")
        painter.setPen(QtGui.QPen(color,3))
        for y in (0.,top):
            painter.drawLine(self._point(camera,(0.,y)),self._point(camera,(right,y)))

    def paintEvent(self, event):
        painter=QtGui.QPainter(self)
        painter.fillRect(self.rect(),QtGui.QColor("#15191f"))
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        camera=self.camera()
        snap=self.snapshot
        if snap.apparatus and self.cylinder_only:
            # The full barrel shows the piston's range of travel.
            painter.setPen(QtGui.QPen(QtGui.QColor("#3a4350"),1,
                                      QtCore.Qt.PenStyle.DashLine))
            width,height=self.barrel
            painter.drawPolyline(QtGui.QPolygonF([self._point(camera,p) for p in
                ((0.,0.),(width,0.),(width,height),(0.,height))]))
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
        if snap.apparatus and self.cylinder_only:
            self._paint_jackets(painter,camera)
        elif snap.apparatus:
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
    save_checkpoint_requested=QtCore.Signal(str)
    load_checkpoint_requested=QtCore.Signal(str)
    export_requested=QtCore.Signal(str)

    def __init__(self,preset,autoplay=True,transient_cycles=2,
                 efficiency_min_cycles=8,shaft_speed=.15,shaft_mode="controlled",
                 config=None):
        super().__init__(); self.setWindowTitle("Microscopic Thermodynamics Laboratory" +
            (" — experimental triangle CCD" if preset=="carnot_triangles" else ""))
        from .lab_view import LabView
        if config is None:
            config=RunConfig(preset=preset,
                             particles=96 if preset=="carnot_triangles" else 48,
                             max_horizon=.02,shaft_speed=shaft_speed,
                             shaft_mode=shaft_mode)
        config.validate()
        preset,shaft_speed,shaft_mode=config.preset,config.shaft_speed,config.shaft_mode
        simulation=load_preset(config)
        self.config=config
        self._generation=0
        self.lab=LabView(simulation)
        self.readout=self.lab.readout
        self.diagrams=self.lab.diagrams
        self.setCentralWidget(self.lab)
        self.lab.setSizes([1100,380])
        self.resize(1480,900)
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
        self.addToolBarBreak()
        bar=self.addToolBar("Run settings")
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
        file_menu=self.menuBar().addMenu("File")
        file_menu.addAction("Save configuration…",lambda: self.save_configuration())
        file_menu.addAction("Open configuration…",lambda: self.open_configuration())
        file_menu.addSeparator()
        file_menu.addAction("Save checkpoint…",lambda: self.save_checkpoint())
        file_menu.addAction("Open checkpoint…",lambda: self.open_checkpoint())
        file_menu.addAction("Export run…",lambda: self.export_run())
        file_menu.addAction("Export plots…",lambda: self.export_plots())
        file_menu.addSeparator()
        self.record_action=file_menu.addAction("Start screen recording…",
                                               lambda: self.start_recording())
        self.stop_record_action=file_menu.addAction("Stop screen recording",
                                                    self.stop_recording)
        self.stop_record_action.setEnabled(False)
        self._record_timer=QtCore.QTimer(self)
        self._record_timer.setInterval(100)
        self._record_timer.timeout.connect(self._capture_recording_frame)
        self._record_folder=None
        self._record_frames=[]
        self._record_start=None
        self._record_writer=None
        self._record_futures=[]
        self._record_dropped=0
        self._record_capture_seconds=[]
        self.status=self.statusBar()
        self.state_label=QtWidgets.QLabel("PAUSED" if not autoplay else "HEALTHY · RUNNING")
        self.status.addPermanentWidget(self.state_label)
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
        self.save_checkpoint_requested.connect(self.worker.save_checkpoint)
        self.load_checkpoint_requested.connect(self.worker.load_checkpoint)
        self.export_requested.connect(self.worker.export)
        self.frame_received.connect(self.worker.frame_received)
        self.stop_requested.connect(self.worker.stop)
        self.worker.snapshot_ready.connect(self.on_snapshot)
        self.worker.failed.connect(self.on_failure)
        self.worker.busy.connect(self.on_busy)
        self.worker.saved.connect(self.on_saved)
        self.worker.checkpoint_loaded.connect(self.on_checkpoint_loaded)
        self.worker.finished.connect(self.thread.quit,
                                     QtCore.Qt.ConnectionType.DirectConnection)
        self.thread.finished.connect(self.worker.deleteLater)
        self.thread.start()

    @property
    def view(self):
        return self.lab.view

    @QtCore.Slot(object)
    def on_snapshot(self, frame):
        if frame.generation != self._generation:
            self.frame_received.emit()
            return
        snapshot=frame.snapshot
        # Paused or stepped frames (no playback rates) always redraw plots.
        series=(InstrumentSeries.from_history(frame.instruments.history)
                if self.lab.plots_due(force=frame.rates is None) else None)
        self.lab.show_frame(snapshot,frame.instruments,series,frame.rates)
        self.status.showMessage(
            f"{'PAUSED' if self.play.isChecked() else 'RUNNING'}  "
            f"t={snapshot.time:.3f}  events={snapshot.event_count}  "
            f"T={snapshot.translational_temperature:.3f}  "
            f"R_E={snapshot.energy_residual:.2e}")
        self.state_label.setText("PAUSED" if self.play.isChecked() else
                                 "HEALTHY · RUNNING")
        self.frame_received.emit()

    @QtCore.Slot(str,str)
    def on_failure(self, message, path):
        self.play.setChecked(True)
        self.state_label.setText("PAUSED")
        detail=f"; diagnostic: {path}" if path else ""
        self.status.showMessage(f"PAUSED: {message}{detail}")
        QtWidgets.QMessageBox.critical(self,"Numerical failure",
            f"The scientific run was paused.\n{message}\n"
            f"Diagnostic and checkpoint: {path}" if path else message)

    @QtCore.Slot()
    def on_busy(self):
        self.state_label.setText("REFINING")
        self.status.showMessage("REFINING  computing the next physics state…")

    @QtCore.Slot(str)
    def on_saved(self, path):
        self.status.showMessage(f"Saved {path}",10000)

    @QtCore.Slot(object,int)
    def on_checkpoint_loaded(self, config, generation):
        self.play.setChecked(True)
        self._install_config(config,generation=generation,reset_worker=False)

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
        self.view.selected_component=None
        self.readout.show_inspection("")
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
        self._install_config(config)

    def _install_config(self, config, generation=None, reset_worker=True):
        self.lab.set_simulation(load_preset(config))
        self.config=config
        self.seed.setValue(config.seed)
        self.shaft_speed.blockSignals(True)
        self.shaft_speed.setValue(config.shaft_speed)
        self.shaft_speed.blockSignals(False)
        self.preset_choice.blockSignals(True)
        self.preset_choice.setCurrentText(config.preset)
        self.preset_choice.blockSignals(False)
        can_drive=config.preset.startswith("carnot_") and config.shaft_mode=="controlled"
        self.branch_action.setEnabled(can_drive)
        self.shaft_speed.setEnabled(can_drive)
        self.setWindowTitle("Microscopic Thermodynamics Laboratory"+
            (" — experimental triangle CCD" if config.preset=="carnot_triangles" else ""))
        self._generation = self._generation+1 if generation is None else generation
        if reset_worker:
            self.configuration_requested.emit(config)

    def save_configuration(self, path=None):
        if path is None:
            path,_=QtWidgets.QFileDialog.getSaveFileName(
                self,"Save configuration","run_config.json","JSON (*.json)")
        if path:
            Path(path).write_text(json.dumps(asdict(self.config),indent=2)+"\n")
            self.on_saved(path)

    def open_configuration(self, path=None):
        if path is None:
            path,_=QtWidgets.QFileDialog.getOpenFileName(
                self,"Open configuration","","JSON (*.json)")
        if not path:
            return
        try:
            data=json.loads(Path(path).read_text())
            names={field.name for field in fields(RunConfig)}
            # Files from older versions may lack fields added since; those
            # take their defaults.
            if not set(data)<=names:
                raise ValueError("configuration fields do not match this version")
            config=RunConfig(**data)
            config.validate()
            self.play.setChecked(True)
            self._install_config(config)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self,"Configuration error",str(exc))

    def save_checkpoint(self, path=None):
        if path is None:
            path,_=QtWidgets.QFileDialog.getSaveFileName(
                self,"Save checkpoint","run_checkpoint.pkl.gz",
                "Compressed checkpoint (*.pkl.gz)")
        if path:
            self.save_checkpoint_requested.emit(str(path))

    def open_checkpoint(self, path=None):
        if path is None:
            path,_=QtWidgets.QFileDialog.getOpenFileName(
                self,"Open locally created checkpoint","",
                "Compressed checkpoint (*.pkl.gz)")
        if path:
            self.play.setChecked(True)
            self.load_checkpoint_requested.emit(str(path))

    def export_run(self, path=None):
        if path is None:
            path=QtWidgets.QFileDialog.getExistingDirectory(self,"Export run to folder")
        if path:
            self.export_requested.emit(str(path))

    def export_plots(self, path=None):
        if path is None:
            path,_=QtWidgets.QFileDialog.getSaveFileName(
                self,"Export plots","plots.png","PNG image (*.png)")
        if path:
            stem=Path(path).with_suffix("")
            for name,plot in (("pressure_area",self.diagrams.pv),
                              ("temperature_entropy",self.diagrams.ts),
                              ("temperatures",self.readout.temperatures)):
                exporter=pyqtgraph.exporters.ImageExporter(plot.plotItem)
                exporter.parameters()["width"]=1200
                exporter.export(str(stem)+f"_{name}.png")
            self.on_saved(str(stem)+"_{pressure_area,temperature_entropy,temperatures}.png")

    def start_recording(self, path=None):
        if path is None:
            path=QtWidgets.QFileDialog.getExistingDirectory(
                self,"Choose folder for screen recording")
        if not path:
            return
        folder=Path(path)
        folder.mkdir(parents=True,exist_ok=True)
        if (folder/"recording.json").exists() or list(folder.glob("frame_*.png")):
            raise ValueError("recording folder already contains frames")
        self._record_folder=folder
        self._record_frames=[]
        self._record_start=time.perf_counter()
        self._record_capture_backend=("xcb-screen" if
            QtGui.QGuiApplication.platformName().lower()=="xcb" else
            "widget")
        self._record_writer=ThreadPoolExecutor(max_workers=1,
                                                thread_name_prefix="screen-record")
        self._record_futures=[]
        self._record_dropped=0
        self._record_capture_seconds=[]
        self.record_action.setEnabled(False)
        self.stop_record_action.setEnabled(True)
        self._capture_recording_frame()
        self._record_timer.start()

    def _capture_recording_frame(self):
        if self._record_folder is None:
            return
        remaining=[]
        for future in self._record_futures:
            if future.done():
                if not future.result():
                    raise OSError("could not save a recording frame")
            else:
                remaining.append(future)
        self._record_futures=remaining
        if len(remaining)>=2:
            self._record_dropped+=1
            return
        capture_start=time.perf_counter()
        if self._record_capture_backend=="xcb-screen":
            screen=self.windowHandle().screen()
            pixmap=screen.grabWindow(int(self.winId())) if screen else QtGui.QPixmap()
            if pixmap.isNull():
                pixmap=self.grab()
                self._record_capture_backend="widget"
        else:
            pixmap=self.grab()
        image=pixmap.toImage()
        self._record_capture_seconds.append(time.perf_counter()-capture_start)
        name=f"frame_{len(self._record_frames):06d}.png"
        self._record_futures.append(self._record_writer.submit(
            image.save,str(self._record_folder/name),"PNG"))
        self._record_frames.append({"file":name,
            "wall_seconds":time.perf_counter()-self._record_start,
            "physical_time":self.view.snapshot.time,
            "event_count":self.view.snapshot.event_count})

    def stop_recording(self):
        if self._record_folder is None:
            return
        self._record_timer.stop()
        for future in self._record_futures:
            if not future.result():
                raise OSError("could not save a recording frame")
        self._record_futures=[]
        if (not self._record_frames or
                self._record_frames[-1]["physical_time"]!=self.view.snapshot.time):
            self._capture_recording_frame()
        folder=self._record_folder
        for future in self._record_futures:
            if not future.result():
                raise OSError("could not save a recording frame")
        self._record_writer.shutdown(wait=True)
        self._record_writer=None
        manifest={"format":"microthermo-screen-frames-v1",
                  "preset":self.config.preset,"seed":self.config.seed,
                  "interval_ms":self._record_timer.interval(),
                  "capture_backend":self._record_capture_backend,
                  "dropped_frames":self._record_dropped,
                  "max_capture_seconds":max(self._record_capture_seconds,
                                            default=0.),
                  "frames":self._record_frames}
        (folder/"recording.json").write_text(json.dumps(manifest,indent=2)+"\n")
        self._record_folder=None
        self.record_action.setEnabled(True)
        self.stop_record_action.setEnabled(False)
        self.on_saved(str(folder/"recording.json"))

    def closeEvent(self, event):
        self.stop_recording()
        self.stop_requested.emit()
        self.thread.wait()
        super().closeEvent(event)


def launch(preset="gas_box",transient_cycles=2,efficiency_min_cycles=8,
           shaft_speed=.15,shaft_mode="controlled",config=None):
    app=QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window=MainWindow(preset,transient_cycles=transient_cycles,
                      efficiency_min_cycles=efficiency_min_cycles,
                      shaft_speed=shaft_speed,shaft_mode=shaft_mode,config=config)
    window.show(); return app.exec()
