"""Qt-thread ownership and bounded snapshot transport for the simulation."""

from __future__ import annotations

from dataclasses import dataclass, replace
from collections import deque
import gzip
import pickle
from pathlib import Path
import tempfile
import time
import traceback
import math
import numpy as np
from PySide6 import QtCore

from ..api import load_preset
from ..runner.edmd_simulation import EdmdCheckpoint
from ..runner.simulation import Checkpoint, NumericalFailure, Snapshot
from ..config import RunConfig
from ..io.exports import export_run
from ..measurements.live import InstrumentFrame, LiveInstruments
from ..measurements.osmosis import OsmosisFrame, OsmosisInstruments


@dataclass(frozen=True)
class RenderFrame:
    snapshot: Snapshot
    instruments: InstrumentFrame
    rates: PlaybackRates | None = None
    generation: int = 0
    osmosis: OsmosisFrame | None = None


@dataclass(frozen=True)
class PlaybackRates:
    physical_per_wall: float
    events_per_wall: float
    wall_window: float


def osmosis_instruments(simulation) -> OsmosisInstruments | None:
    """Osmosis instruments for worlds with ring hosts behind a membrane."""
    metadata = simulation.world.metadata
    return OsmosisInstruments() if "membrane_x" in metadata and simulation.world.rings else None


def immutable_snapshot(snapshot: Snapshot) -> Snapshot:
    """Detach arrays into bytes-backed, genuinely read-only render buffers."""
    def frozen(array):
        return np.frombuffer(array.tobytes(),dtype=array.dtype).reshape(array.shape)
    return replace(snapshot, position=frozen(snapshot.position),
                   velocity=frozen(snapshot.velocity),angle=frozen(snapshot.angle),
                   omega=frozen(snapshot.omega),
                   membership=(None if snapshot.membership is None
                               else frozen(snapshot.membership)))


class SimulationWorker(QtCore.QObject):
    snapshot_ready = QtCore.Signal(object)
    failed = QtCore.Signal(str, str)
    finished = QtCore.Signal()
    busy = QtCore.Signal()
    saved = QtCore.Signal(str)
    checkpoint_loaded = QtCore.Signal(object, int)

    def __init__(self, simulation, autoplay=True, physical_step=.016,
                 transient_cycles=2, efficiency_min_cycles=8, config=None):
        super().__init__()
        self._simulation=simulation
        self._autoplay=autoplay
        self._physical_step=physical_step
        self._base_step=physical_step
        self._config=config
        self._generation=0
        self._transient_cycles=transient_cycles
        self._efficiency_min_cycles=efficiency_min_cycles
        self._timer=None
        self._osmosis=osmosis_instruments(simulation)
        self._instruments=LiveInstruments(
            history_limit=3000,history_spacing=.05,
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
        osmosis=(self._osmosis.observe(self._simulation,snapshot)
                 if self._osmosis is not None else None)
        self.snapshot_ready.emit(RenderFrame(immutable_snapshot(snapshot),instruments,
                                             self._rates,self._generation,osmosis))

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
            self.busy.emit()
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

    @QtCore.Slot(str)
    def save_checkpoint(self, path):
        try:
            payload={"config":self._config,
                     "checkpoint":self._simulation.checkpoint()}
            target=Path(path)
            with gzip.open(target,"wb",compresslevel=1) as stream:
                pickle.dump(payload,stream,protocol=pickle.HIGHEST_PROTOCOL)
            self.saved.emit(str(target))
        except Exception as exc:
            self.failed.emit(f"Checkpoint save failed: {exc}","")

    @QtCore.Slot(str)
    def load_checkpoint(self, path):
        try:
            # Pickle is restricted to files deliberately selected by the user.
            with gzip.open(path,"rb") as stream:
                payload=pickle.load(stream)
            config=payload["config"]
            checkpoint=payload["checkpoint"]
            if not isinstance(config,RunConfig) or not isinstance(checkpoint,(Checkpoint,EdmdCheckpoint)):
                raise ValueError("unrecognized checkpoint")
            simulation=load_preset(config)
            simulation.restore(checkpoint)
            if self._timer is not None:
                self._timer.stop()
            self._simulation=simulation
            self._config=config
            self._rates=None
            self._osmosis=osmosis_instruments(simulation)
            self._instruments=LiveInstruments(
                history_limit=3000,history_spacing=.05,
                transient_cycles=self._transient_cycles,
                efficiency_min_cycles=self._efficiency_min_cycles)
            self._generation += 1
            self._awaiting_frame=False
            self._dirty_frame=False
            self.checkpoint_loaded.emit(config,self._generation)
            self._publish()
        except Exception as exc:
            self.failed.emit(f"Checkpoint load failed: {exc}","")

    @QtCore.Slot(str)
    def export(self, path):
        try:
            export_run(path,self._config,self._simulation)
            self.saved.emit(str(path))
        except Exception as exc:
            self.failed.emit(f"Export failed: {exc}","")

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

    @QtCore.Slot(float)
    def set_playback_rate(self, rate):
        if not math.isfinite(rate) or rate <= 0:
            raise ValueError("playback rate must be positive and finite")
        self._physical_step=self._base_step*rate
        if self._timer is not None and self._timer.isActive():
            self._reset_rates()

    @QtCore.Slot()
    def step_branch(self):
        mechanism=self._simulation.world.mechanism
        shaft=getattr(mechanism,"shaft",None)
        cam=getattr(mechanism,"cam",None)
        if shaft is None or cam is None or shaft.prescribed_omega is None:
            return
        if self._timer is not None:
            self._timer.stop()
        self._rates=None
        sign=-1. if mechanism.reversed_cycle else 1.
        phase=sign*shaft.phi
        rate=sign*shaft.prescribed_omega
        if rate <= 0:
            return
        period=2*math.pi
        cycle=math.floor(phase/period)
        candidates=[cycle*period+float(boundary) for boundary in
                    cam.boundaries[1:]]
        candidates.append((cycle+1)*period+float(cam.boundaries[1]))
        target_phase=next(value for value in candidates if value>phase+1e-10)
        target_time=self._simulation.time+(target_phase-phase)/rate
        self._run(lambda: self._simulation.advance_to(target_time))

    @QtCore.Slot(float)
    def set_shaft_speed(self, speed):
        if self._timer is not None:
            self._timer.stop()
        self._rates=None
        def change():
            self._simulation.apply_command({"name":"set_shaft_speed",
                                            "speed":speed})
            if self._config is not None:
                self._config=replace(self._config,shaft_speed=speed)
        self._run(change)

    @QtCore.Slot(float)
    def set_wall_temperature(self, temperature):
        """Change the thermal walls' temperature now (an intervention)."""
        def change():
            self._simulation.apply_command({"name":"set_wall_temperature",
                                            "temperature":temperature})
        self._run(change)

    @QtCore.Slot(int)
    def reset_seed(self, seed):
        if self._config is None:
            return
        self._reset_configuration(replace(self._config,seed=seed))

    @QtCore.Slot(object)
    def reset_configuration(self, config):
        self._reset_configuration(config)

    def _reset_configuration(self, config):
        if self._stopping:
            return
        if self._timer is not None:
            self._timer.stop()
        self._rates=None
        try:
            simulation=load_preset(config)
            self._config=config
            self._simulation=simulation
            self._osmosis=osmosis_instruments(simulation)
            self._instruments=LiveInstruments(
                history_limit=3000,history_spacing=.05,
                transient_cycles=self._transient_cycles,
                efficiency_min_cycles=self._efficiency_min_cycles)
            self._generation += 1
            self._awaiting_frame=False
            self._dirty_frame=False
            self._publish()
        except Exception as exc:
            self.failed.emit(f"{type(exc).__name__}: {exc}","")

    @QtCore.Slot()
    def stop(self):
        if self._stopping:
            return
        self._stopping=True
        if self._timer is not None:
            self._timer.stop()
        self.finished.emit()
