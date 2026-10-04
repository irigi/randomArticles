"""Time-based desktop playback of a precalculated archive."""

from __future__ import annotations

import gc
import time

from PySide6 import QtCore, QtWidgets

from ..api import load_preset
from ..config import RunConfig
from ..io.replay import ReplayReader
from .main_window import ApparatusView


class ReplayWindow(QtWidgets.QMainWindow):
    def __init__(self, path):
        super().__init__()
        self.reader = ReplayReader(path)
        self.reader.prepare_playback()
        config = RunConfig(**self.reader.manifest["config"])
        self.view = ApparatusView(load_preset(config))
        self.duration = self.reader.manifest["duration"]
        self.physical_time = 0.0
        self.playback_speed = 1.0
        self._clock = None
        self._seeking = False
        self._restore_gc = False
        self.setWindowTitle(f"Replay — {config.preset} — seed {config.seed}")
        root = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(root)
        layout.addWidget(self.view, stretch=1)
        controls = QtWidgets.QHBoxLayout()
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
        layout.addLayout(controls)
        self.readout = QtWidgets.QLabel()
        self.readout.setWordWrap(True)
        layout.addWidget(self.readout)
        self.setCentralWidget(root)
        self.resize(1100, 700)
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

    def _render(self):
        snapshot, index = self.reader.sample(self.physical_time)
        self.view.set_snapshot(snapshot)
        ledger = self.reader.ledger(index)
        recorded = self.reader.frame(index)
        self.readout.setText(
            f"t={snapshot.time:.3f}/{self.duration:.3f} s · "
            f"recorded statistics at t={recorded.time:.3f} s · "
            f"events={recorded.event_count} · "
            f"T={recorded.translational_temperature:.3f} · "
            f"QH={ledger['heat_hot']:.3f} · QC={ledger['heat_cold']:.3f} · "
            f"motor work={ledger['work_on']:.3f} · "
            f"load={ledger['load_output']:.3f} · "
            f"R_E={recorded.energy_residual:.2e}")
        self._seeking = True
        self.slider.setValue(round(10000*self.physical_time/self.duration))
        self._seeking = False


def launch_replay(path):
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = ReplayWindow(path)
    window.show()
    return app.exec()
