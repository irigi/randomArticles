"""Qt-thread ownership and bounded snapshot transport for the simulation."""

from __future__ import annotations

from dataclasses import dataclass, replace
from collections import deque
import tempfile
import time
import traceback
import numpy as np
from PySide6 import QtCore

from ..runner.simulation import NumericalFailure, Snapshot
from ..measurements.live import InstrumentFrame, LiveInstruments


@dataclass(frozen=True)
class RenderFrame:
    snapshot: Snapshot
    instruments: InstrumentFrame
    rates: PlaybackRates | None = None


@dataclass(frozen=True)
class PlaybackRates:
    physical_per_wall: float
    events_per_wall: float
    wall_window: float


def immutable_snapshot(snapshot: Snapshot) -> Snapshot:
    """Detach arrays into bytes-backed, genuinely read-only render buffers."""
    def frozen(array):
        return np.frombuffer(array.tobytes(),dtype=array.dtype).reshape(array.shape)
    return replace(snapshot, position=frozen(snapshot.position),
                   velocity=frozen(snapshot.velocity),angle=frozen(snapshot.angle),
                   omega=frozen(snapshot.omega))


class SimulationWorker(QtCore.QObject):
    snapshot_ready = QtCore.Signal(object)
    failed = QtCore.Signal(str, str)
    finished = QtCore.Signal()

    def __init__(self, simulation, autoplay=True, physical_step=.01,
                 transient_cycles=2, efficiency_min_cycles=8):
        super().__init__()
        self._simulation=simulation
        self._autoplay=autoplay
        self._physical_step=physical_step
        self._timer=None
        self._instruments=LiveInstruments(
            transient_cycles=transient_cycles,
            efficiency_min_cycles=efficiency_min_cycles)
        self._awaiting_frame=False
        self._dirty_frame=False
        self._stopping=False
        self._rate_clock=None
        self._rate_samples=deque()
        self._rates=None

    def _reset_rates(self):
        self._rate_samples.clear()
        self._rates=None
        self._rate_clock=time.perf_counter()

    def _record_rate(self, physical_delta, event_delta):
        now=time.perf_counter()
        elapsed=max(now-self._rate_clock,1e-12)
        self._rate_clock=now
        self._rate_samples.append((elapsed,physical_delta,event_delta))
        while len(self._rate_samples)>1 and sum(s[0] for s in self._rate_samples)>1.:
            self._rate_samples.popleft()
        wall=sum(s[0] for s in self._rate_samples)
        self._rates=PlaybackRates(
            sum(s[1] for s in self._rate_samples)/wall,
            sum(s[2] for s in self._rate_samples)/wall,wall)

    @QtCore.Slot()
    def start(self):
        if self._stopping:
            return
        self._timer=QtCore.QTimer(self)
        self._timer.setInterval(16)
        self._timer.timeout.connect(self.tick)
        self._publish()
        if self._autoplay:
            self._reset_rates()
            self._timer.start()

    def _publish(self):
        if self._awaiting_frame:
            self._dirty_frame=True
            return
        self._awaiting_frame=True
        self._dirty_frame=False
        snapshot=self._simulation.snapshot()
        instruments=self._instruments.observe(self._simulation,snapshot)
        self.snapshot_ready.emit(RenderFrame(immutable_snapshot(snapshot),instruments,
                                             self._rates))

    @QtCore.Slot()
    def frame_received(self):
        self._awaiting_frame=False
        if self._dirty_frame and not self._stopping:
            self._publish()

    @QtCore.Slot(bool)
    def set_playing(self, playing):
        if self._stopping or self._timer is None:
            return
        if playing:
            self._reset_rates()
            self._timer.start()
        else:
            self._timer.stop()
            self._rates=None
            self._publish()

    def _run(self, action, timed=False):
        if self._stopping:
            return
        try:
            before_time=self._simulation.time
            before_events=self._simulation.event_count
            action()
            if timed:
                self._record_rate(self._simulation.time-before_time,
                                  self._simulation.event_count-before_events)
            self._publish()
        except NumericalFailure as exc:
            self._timer.stop()
            path=self._simulation.save_failure(
                tempfile.mkdtemp(prefix="microthermo-failure-"))
            self.failed.emit(str(exc),str(path))
        except Exception as exc:
            self._timer.stop()
            self.failed.emit(f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}","")

    @QtCore.Slot()
    def tick(self):
        self._run(lambda: self._simulation.advance_to(
            self._simulation.time+self._physical_step),timed=True)

    @QtCore.Slot()
    def step_collision(self):
        if self._timer is not None:
            self._timer.stop()
        self._rates=None
        self._run(self._simulation.step_collision)

    @QtCore.Slot(float)
    def advance_duration(self, duration):
        if duration > 0:
            self._rates=None
            self._run(lambda: self._simulation.advance_to(
                self._simulation.time+duration))

    @QtCore.Slot()
    def stop(self):
        if self._stopping:
            return
        self._stopping=True
        if self._timer is not None:
            self._timer.stop()
        self.finished.emit()
