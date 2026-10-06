"""Time-based desktop playback of a precalculated archive."""

from __future__ import annotations

import gc
import time

from PySide6 import QtCore, QtWidgets

from ..api import load_preset
from ..config import RunConfig
from ..io.replay import ReplayReader
from ..measurements.ideal import ideal_carnot_reference
from ..measurements.live import InstrumentFrame, gas_degrees_of_freedom
from .instruments import InstrumentSeries
from .lab_view import LabView


class ReplayWindow(QtWidgets.QMainWindow):
    def __init__(self, path):
        super().__init__()
        self.reader = ReplayReader(path)
        self.reader.prepare_playback()
        config = RunConfig(**self.reader.manifest["config"])
        simulation = load_preset(config)
        self.lab = LabView(simulation)
        self.instruments = self.reader.has_instruments
        if self.instruments:
            self._columns = self.reader.instrument_series()
            self._reference = ideal_carnot_reference(simulation)
            self._particles = simulation.world.bodies.n
            self._dof = gas_degrees_of_freedom(simulation)
        else:
            # Version-1 archives hold no instrument history; show the cylinder
            # and the one-line ledger only.
            self.lab.diagrams.setVisible(False)
            self.lab.readout.setVisible(False)
        self.osmosis = self.reader.has_osmosis
        self._osmosis_frame = None
        self._shown_index = None
        self._plotted_index = None
        self.duration = self.reader.manifest["duration"]
        self.physical_time = 0.0
        self.playback_speed = 1.0
        self._clock = None
        self._seeking = False
        self._restore_gc = False
        self.setWindowTitle(f"Replay — {config.preset} — seed {config.seed}")
        root = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(root)
        layout.addWidget(self.lab, stretch=1)
        self.controls = QtWidgets.QWidget()
        controls = QtWidgets.QHBoxLayout(self.controls)
        controls.setContentsMargins(0, 0, 0, 0)
        self.play = QtWidgets.QPushButton("Play")
        self.play.setCheckable(True)
        self.play.toggled.connect(self._set_playing)
        controls.addWidget(self.play)
        self.speed = QtWidgets.QComboBox()
        for rate in (.25, .5, 1., 2., 4.):
            self.speed.addItem(f"{rate:g}×", rate)
        self.speed.setCurrentIndex(2)
        self.speed.currentIndexChanged.connect(self._set_speed)
        controls.addWidget(self.speed)
        self.slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.slider.setRange(0, 10000)
        self.slider.sliderPressed.connect(self._pause_for_seek)
        self.slider.valueChanged.connect(self._seek_slider)
        controls.addWidget(self.slider, stretch=1)
        layout.addWidget(self.controls)
        self.status_line = QtWidgets.QLabel()
        self.status_line.setWordWrap(True)
        layout.addWidget(self.status_line)
        self.setCentralWidget(root)
        self.lab.setSizes([1100, 380])
        self.resize(1480, 940)
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(16)
        self.timer.timeout.connect(self._tick)
        self._render()

    def closeEvent(self, event):
        self.timer.stop()
        self._restore_gc_state()
        self.reader.close()
        super().closeEvent(event)

    def _restore_gc_state(self):
        if self._restore_gc:
            gc.enable()
            self._restore_gc = False

    def _set_playing(self, playing):
        self.play.setText("Pause" if playing else "Play")
        if playing:
            if self.physical_time >= self.duration:
                self.physical_time = 0.
            if gc.isenabled():
                gc.disable()
                self._restore_gc = True
            self._clock = time.perf_counter()
            self.timer.start()
        else:
            self.timer.stop()
            self._clock = None
            self._restore_gc_state()
            self._render()

    def _set_speed(self):
        self.playback_speed = float(self.speed.currentData())
        self._clock = time.perf_counter() if self.play.isChecked() else None

    def _pause_for_seek(self):
        self.play.setChecked(False)

    def _seek_slider(self, value):
        if self._seeking:
            return
        self.physical_time = self.duration * value / 10000
        self._render()

    def _tick(self):
        now = time.perf_counter()
        self.physical_time = min(
            self.duration, self.physical_time +
            (now-self._clock)*self.playback_speed)
        self._clock = now
        self._render()
        if self.physical_time >= self.duration:
            self.play.setChecked(False)

    @property
    def view(self):
        return self.lab.view

    def _instrument_frame(self, index):
        current = self.reader.instrument_sample(index)
        summaries, _, markers = self.reader.cycles()
        completed = current.completed_cycles
        return InstrumentFrame(
            current, (), self._reference, summaries[:completed],
            markers[completed] if completed < len(markers) else None,
            self._particles, self._dof)

    def _render(self):
        snapshot, index = self.reader.sample(self.physical_time)
        if self.instruments:
            if index != self._shown_index:
                self._shown_index = index
                self._frame = self._instrument_frame(index)
                if self.osmosis:
                    self._osmosis_frame = self.reader.osmosis_frame(index)
                self._plotted_index = None
            series = None
            if (self._plotted_index != index and
                    self.lab.plots_due(force=not self.play.isChecked())):
                self._plotted_index = index
                series = InstrumentSeries.from_columns(self._columns, index+1,
                                                       max_points=3000)
            self.lab.show_frame(snapshot, self._frame, series, replay=True,
                                osmosis=self._osmosis_frame)
        else:
            self.lab.view.set_snapshot(snapshot)
            if self.lab.dial is not None:
                self.lab.dial.set_snapshot(snapshot)
        ledger = self.reader.ledger(index)
        recorded = self.reader.frame(index)
        if self._osmosis_frame is not None:
            energy = (f"wall T={self._osmosis_frame.current.temperature:.3f} · "
                      f"heat in from the walls={ledger['heat_hot']:.3f} · ")
        else:
            energy = (f"QH={ledger['heat_hot']:.3f} · QC={ledger['heat_cold']:.3f} · "
                      f"motor work={ledger['work_on']:.3f} · "
                      f"load={ledger['load_output']:.3f} · ")
        self.status_line.setText(
            f"t={snapshot.time:.3f}/{self.duration:.3f} s · "
            f"recorded statistics at t={recorded.time:.3f} s · "
            f"events={recorded.event_count} · "
            f"T={recorded.translational_temperature:.3f} · {energy}"
            f"R_E={recorded.energy_residual:.2e}")
        self._seeking = True
        self.slider.setValue(round(10000*self.physical_time/self.duration))
        self._seeking = False


def launch_replay(path):
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = ReplayWindow(path)
    window.show()
    return app.exec()
