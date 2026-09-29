from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets
import math
from ..api import load_preset
from ..config import RunConfig
from .instruments import InstrumentPanel
from .scene import Camera, scene_bounds
from .worker import SimulationWorker


class ApparatusView(QtWidgets.QWidget):
    def __init__(self, sim):
        super().__init__()
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
        self.setMinimumSize(500,340)

    def set_snapshot(self, snapshot):
        self.snapshot=snapshot
        self.update()

    def _point(self, camera, point):
        return QtCore.QPointF(*camera.map(*point))

    def _component(self, painter, camera, component):
        colors = {
            "cylinder":"#8997a7", "piston":"#f0c46c",
            "piston_cam":"#75bdd5", "cam_follower":"#e8d8a8",
            "selector_cam":"#ab9edb", "selector_shoe":"#e8d8a8",
            "thermal_wall":{"hot":"#f2a65a", "cold":"#76bce3"}.get(
                component.state,"#d5dce4"),
            "flywheel":"#d8dee7", "flywheel_spoke":"#f0c46c",
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
        camera=Camera(self.bounds,self.width(),self.height())
        snap=self.snapshot
        if snap.apparatus:
            for component in snap.apparatus:
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
        if snap.apparatus:
            painter.setPen(QtGui.QColor("#d5dce4"))
            parts={part.name:part for part in snap.apparatus}
            def center(name):
                points=parts[name].points
                return ((min(x for x,_ in points)+max(x for x,_ in points))/2,
                        (min(y for _,y in points)+max(y for _,y in points))/2)
            shaft_center=parts["flywheel_spoke"].points[0]
            labels=(("Piston",parts["piston"].points[1],-50,-8),
                    ("Flywheel",shaft_center,-26,43),
                    ("Piston cam",center("piston_cam"),-36,43),
                    ("Selector cam",center("selector_cam"),-40,43),
                    ("Spring",shaft_center,-19,-33),
                    ("Load",parts["load"].points[-1],7,5))
            for label,point,offset_x,offset_y in labels:
                screen=self._point(camera,point)
                painter.drawText(screen+QtCore.QPointF(offset_x,offset_y),label)
            painter.drawText(16,24,f"{snap.branch.replace('_',' ')} · phase {snap.shaft_phase:.2f} rad")
        painter.end()


class MainWindow(QtWidgets.QMainWindow):
    play_requested=QtCore.Signal(bool)
    step_requested=QtCore.Signal()
    advance_requested=QtCore.Signal(float)
    frame_received=QtCore.Signal()
    stop_requested=QtCore.Signal()

    def __init__(self,preset,autoplay=True,transient_cycles=2,
                 efficiency_min_cycles=8,shaft_speed=.15):
        super().__init__(); self.setWindowTitle("Microscopic Thermodynamics Laboratory" +
            (" — experimental triangle CCD" if preset=="carnot_triangles" else ""))
        simulation=load_preset(RunConfig(preset=preset,particles=48,
                                         max_horizon=.02,shaft_speed=shaft_speed))
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
        self.play=bar.addAction("Pause"); self.play.setCheckable(True)
        self.play.setChecked(not autoplay)
        self.play.toggled.connect(lambda paused: self.play_requested.emit(not paused))
        step=bar.addAction("Step collision"); step.triggered.connect(self.step_collision)
        self.status=self.statusBar()
        self.thread=QtCore.QThread(self)
        self.worker=SimulationWorker(
            simulation,autoplay=autoplay,transient_cycles=transient_cycles,
            efficiency_min_cycles=efficiency_min_cycles)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.start)
        self.play_requested.connect(self.worker.set_playing)
        self.step_requested.connect(self.worker.step_collision)
        self.advance_requested.connect(self.worker.advance_duration)
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
        snapshot=frame.snapshot
        self.view.set_snapshot(snapshot)
        self.instruments.set_frame(frame.instruments)
        self.instruments.set_rates(frame.rates)
        self.status.showMessage(
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

    def closeEvent(self, event):
        self.stop_requested.emit()
        self.thread.wait()
        super().closeEvent(event)


def launch(preset="gas_box",transient_cycles=2,efficiency_min_cycles=8,
           shaft_speed=.15):
    app=QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window=MainWindow(preset,transient_cycles=transient_cycles,
                      efficiency_min_cycles=efficiency_min_cycles,
                      shaft_speed=shaft_speed)
    window.show(); return app.exec()
