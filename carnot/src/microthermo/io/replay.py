"""Lossless, chunked snapshot archive for offline playback (Carnot and osmosis)."""

from __future__ import annotations

from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, replace
import bisect
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from ..api import load_preset
from ..config import RunConfig
from ..core.apparatus import ApparatusComponent
from ..measurements.cycles import CycleMarker, assess_efficiency, summarize_cycle
from ..measurements.live import InstrumentSample, LiveInstruments
from ..measurements.osmosis import OsmosisRecorder, RecordedOsmosis
from ..measurements.pressure_area import PressureAreaComparison
from ..runner.simulation import Snapshot
from .exports import _json


FORMAT_VERSION = 3
SUPPORTED_VERSIONS = (1, 2, 3)
# Version 2 adds per-frame instrument columns (prefixed "i_") and the cycle
# markers with their pressure-area comparisons in the manifest, so a replay
# can show the same instruments as the live window. Version 3 adds, for
# osmosis worlds, the well membership of every body per frame and the
# osmosis instrument columns (prefixed "o_", see measurements.osmosis);
# Carnot archives only change their version number.
INSTRUMENT_COLUMNS = {
    "area": np.float64, "pressure": np.float64, "gas_energy": np.float64,
    "entropy": np.float64, "motor_work": np.float64,
    "piston_energy": np.float64, "flywheel_energy": np.float64,
    "spring_energy": np.float64, "completed_cycles": np.int64,
    "max_penetration": np.float64, "ccd_refinements": np.int64,
    "ccd_failures": np.int64, "cluster_count": np.int64,
    "max_cluster_residual": np.float64,
}


def _interpolate_apparatus(left: Snapshot, right: Snapshot, alpha: float,
                           branch_switched: bool) -> tuple[ApparatusComponent, ...]:
    apparatus = []
    for a, b in zip(left.apparatus, right.apparatus):
        if a.name != b.name or len(a.points) != len(b.points):
            apparatus.append(a)
            continue
        if a.name == "selector_shoe" and a.state != b.state:
            # The thermal selector jumps between sector detents. Its position
            # changes at the logged branch event, not gradually across a frame.
            points = b.points if branch_switched else a.points
        else:
            points = tuple((x[0]+(y[0]-x[0])*alpha,
                            x[1]+(y[1]-x[1])*alpha)
                           for x, y in zip(a.points, b.points))
        apparatus.append(replace(a, points=points,
                                 state=b.state if branch_switched else a.state))
    return tuple(apparatus)


def _frame_arrays(frames: list[Snapshot]) -> dict[str, np.ndarray]:
    """Keep all physical coordinates and observables in float64."""
    columns = {
        "time": np.array([s.time for s in frames], dtype=np.float64),
        "position": np.stack([s.position for s in frames]),
        "velocity": np.stack([s.velocity for s in frames]),
        "angle": np.stack([s.angle for s in frames]),
        "omega": np.stack([s.omega for s in frames]),
        "energy": np.array([s.energy for s in frames], dtype=np.float64),
        "temperature_trans": np.array(
            [s.translational_temperature for s in frames], dtype=np.float64),
        "temperature_rot": np.array(
            [s.rotational_temperature for s in frames], dtype=np.float64),
        "event_count": np.array([s.event_count for s in frames], dtype=np.int64),
        "energy_residual": np.array(
            [s.energy_residual for s in frames], dtype=np.float64),
        "shaft_phase": np.array(
            [math.nan if s.shaft_phase is None else s.shaft_phase for s in frames],
            dtype=np.float64),
        "occupancy": np.array([s.occupancy for s in frames], dtype=np.int64),
        "branch": np.array([s.branch or "" for s in frames], dtype="U32"),
        "apparatus": np.array([json.dumps([asdict(x) for x in s.apparatus])
                               for s in frames], dtype="U"),
    }
    if frames[0].membership is not None:
        columns["membership"] = np.stack([s.membership for s in frames]).astype(np.int32)
    return columns


def _is_osmosis(simulation) -> bool:
    return bool(simulation.world.rings) and "membrane_x" in simulation.world.metadata


def precalculate_replay(config, output: str | Path, *, fps: float = 60.,
                        chunk_frames: int = 256, progress=None) -> Path:
    """Simulate once and flush bounded frame chunks; store exact event records.

    ``progress(physical_time, duration, event_count)`` is called after each
    chunk is written.
    """
    if not math.isfinite(fps) or fps <= 0 or chunk_frames < 2:
        raise ValueError("fps must be positive and chunk_frames at least two")
    config.validate()
    if not math.isfinite(config.duration):
        raise ValueError("replay duration must be finite")
    path = Path(output)
    path.mkdir(parents=True, exist_ok=True)
    if any(path.iterdir()):
        raise FileExistsError(f"replay output directory is not empty: {path}")
    simulation = load_preset(config)
    bodies = simulation.world.bodies
    geometry = {
        "shape": bodies.shape.tolist(),
        "mass": bodies.mass.tolist(),
        "inertia": bodies.inertia.tolist(),
        "dynamic": bodies.dynamic.tolist(),
        "radius": bodies.radius.tolist(),
        "polygons": {str(k): v.tolist() for k, v in bodies.polygons.items()},
        "walls": [{"name": wall.name, "point": wall.point_at(0.).tolist(),
                   "normal": wall.inward_normal.tolist(),
                   "length": wall.length, "initial_kind": wall.kind_at(0.).value}
                  for wall in simulation.world.walls],
    }
    frames: list[Snapshot] = []
    ledgers: list[tuple[float, float, float, float]] = []
    samples: list[InstrumentSample] = []
    instruments = LiveInstruments(history_limit=2,
                                  transient_cycles=config.transient_cycles,
                                  efficiency_min_cycles=config.efficiency_min_cycles)
    osmosis = OsmosisRecorder() if _is_osmosis(simulation) else None
    chunks = []
    frame_count = 0

    def flush() -> None:
        nonlocal frame_count
        if not frames:
            return
        name = f"frames_{len(chunks):06d}.npz"
        columns = _frame_arrays(frames)
        columns["ledger"] = np.array(ledgers, dtype=np.float64)
        for column, dtype in INSTRUMENT_COLUMNS.items():
            values = [getattr(sample, column) for sample in samples]
            columns["i_"+column] = np.array(
                [math.nan if value is None else value for value in values], dtype=dtype)
        if osmosis is not None:
            columns.update(osmosis.take_columns())
        np.savez_compressed(path/name, **columns)
        chunks.append({"file": name, "first_frame": frame_count,
                       "count": len(frames), "first_time": frames[0].time,
                       "last_time": frames[-1].time})
        frame_count += len(frames)
        if progress is not None:
            progress(frames[-1].time, config.duration, frames[-1].event_count)
        frames.clear()
        ledgers.clear()
        samples.clear()

    def append_frame(snapshot: Snapshot, events) -> None:
        frames.append(snapshot)
        samples.append(instruments.observe(simulation, snapshot, events=events).current)
        if osmosis is not None:
            osmosis.observe(simulation, snapshot)
        ledger = simulation.ledger
        ledgers.append((ledger.heat_hot.value, ledger.heat_cold.value,
                        ledger.work_on.value, ledger.load_output.value))

    append_frame(simulation.snapshot(), simulation.drain_events())
    sample = 1
    max_kind_length, max_participant_length = 1, 2
    with (path/"events.jsonl").open("w") as file:
        while True:
            target = min(config.duration, sample/fps)
            if target <= simulation.time + 1e-12:
                break
            snapshot = simulation.advance_to(target)
            events = simulation.drain_events()
            append_frame(snapshot, events)
            for event in events:
                file.write(json.dumps(event.as_dict(), default=_json)+"\n")
                max_kind_length = max(max_kind_length, len(event.kind))
                max_participant_length = max(
                    max_participant_length, len(json.dumps(event.participants)))
            if len(frames) >= chunk_frames:
                flush()
            if target >= config.duration - 1e-12:
                break
            sample += 1
    flush()
    # Build the original v1 index with disk-backed columns. This keeps the
    # export's memory independent of the number of collision records.
    count = simulation.event_count
    with TemporaryDirectory(dir=path) as scratch:
        def column(name, dtype, shape):
            if count == 0:
                return np.empty(shape, dtype=dtype)
            return np.lib.format.open_memmap(Path(scratch)/f"{name}.npy",
                                             mode="w+", dtype=dtype, shape=shape)

        times = column("time", np.float64, (count,))
        impulses = column("impulse", np.float64, (count, 2))
        points = column("point", np.float64, (count, 2))
        kinds = column("kind", f"U{max_kind_length}", (count,))
        participants = column("participants", f"U{max_participant_length}",
                              (count,))
        with (path/"events.jsonl").open() as file:
            written = 0
            for written, line in enumerate(file, start=1):
                event = json.loads(line)
                index = written-1
                times[index] = event["time"]
                impulses[index] = event["impulse"]
                points[index] = event.get("contact_point") or (math.nan, math.nan)
                kinds[index] = event["kind"]
                participants[index] = json.dumps(event["participants"])
        if written != count:
            raise ValueError("replay event log and event count disagree")
        np.savez_compressed(path/"event_times.npz", time=times, impulse=impulses,
                            point=points, kind=kinds, participants=participants)
    final = simulation.snapshot()
    manifest = {
        "format": "microthermo-replay", "version": FORMAT_VERSION,
        "config": asdict(config), "config_hash": config.digest(),
        "seed": config.seed, "backend": simulation.numeric_backend,
        "sample_fps": fps, "duration": final.time,
        "frame_count": frame_count, "event_count": simulation.event_count,
        "event_times_file": "event_times.npz",
        "chunk_frames": chunk_frames, "chunks": chunks, "geometry": geometry,
        "world_metadata": {k: v for k, v in simulation.world.metadata.items()
                           if k != "clock"},
        "final_ledger": asdict(simulation.ledger),
        "final_energy_residual": final.energy_residual,
        "reconstruction": "float64 snapshots; exact at saved timestamps",
        "ledger_columns": ["heat_hot", "heat_cold", "work_on", "load_output"],
        "instruments": {
            "pressure_window": instruments.pressure_window,
            "transient_cycles": config.transient_cycles,
            "efficiency_min_cycles": config.efficiency_min_cycles,
            "cycle_markers": [asdict(m) for m in simulation.cycle_markers],
            "cycle_pressure_area": [None if c is None else asdict(c)
                                    for c in instruments.cycle_pressure_area],
        },
    }
    if osmosis is not None:
        manifest["osmosis"] = osmosis.settings(simulation)
    (path/"manifest.json").write_text(json.dumps(manifest, indent=2,
                                                 default=_json)+"\n")
    return path


class ReplayReader:
    """Bounded-memory random access to saved snapshots."""

    CACHE_CHUNKS = 3

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.manifest = json.loads((self.path/"manifest.json").read_text())
        if (self.manifest.get("format") != "microthermo-replay" or
                self.manifest.get("version") not in SUPPORTED_VERSIONS):
            raise ValueError("unsupported replay format or version")
        self._series = None
        self._osmosis = None
        self._cycles = None
        self._efficiency = {}
        self.chunks = self.manifest["chunks"]
        self._starts = [chunk["first_frame"] for chunk in self.chunks]
        self._data = None
        self._chunk_index = None
        self._cache = OrderedDict()
        self._executor = None
        self._pending_index = None
        self._pending = None
        self._event_times = None
        self._event_kinds = None
        self._event_participants = None
        self._event_a = None
        self._event_b = None
        self._event_impulse = None
        self._event_point = None
        self._mass = None
        self._inertia = None
        self._dynamic = None
        self._exact_disc = all(shape == 0 for shape in
                               self.manifest["geometry"]["shape"])
        self._exact_motion = None

    def __len__(self) -> int:
        return self.manifest["frame_count"]

    def _read_chunk(self, chunk_index: int) -> dict[str, np.ndarray]:
        with np.load(self.path/self.chunks[chunk_index]["file"],
                     allow_pickle=False) as archive:
            return {name: archive[name] for name in archive.files}

    def _remember_chunk(self, chunk_index: int, data: dict[str, np.ndarray]) -> None:
        self._cache[chunk_index] = data
        self._cache.move_to_end(chunk_index)
        while len(self._cache) > self.CACHE_CHUNKS:
            self._cache.popitem(last=False)

    def _prefetch_next(self, chunk_index: int) -> None:
        next_index = chunk_index+1
        if (self._executor is None or next_index >= len(self.chunks) or
            next_index in self._cache):
            return
        if self._pending is not None:
            if self._pending_index == next_index:
                return
            if not self._pending.done():
                return
            self._pending = None
            self._pending_index = None
        self._pending_index = next_index
        self._pending = self._executor.submit(self._read_chunk, next_index)

    def prepare_playback(self) -> None:
        """Pay initial I/O cost before showing the window; prefetch later chunks."""
        self._times()
        if self._can_reconstruct_motion():
            self._body_properties()
        for chunk_index in range(min(self.CACHE_CHUNKS, len(self.chunks))):
            if chunk_index not in self._cache:
                self._remember_chunk(chunk_index, self._read_chunk(chunk_index))
        if self._executor is None:
            self._executor = ThreadPoolExecutor(max_workers=1,
                                                thread_name_prefix="replay-prefetch")
        self._load(0)

    def close(self) -> None:
        if self._executor is not None:
            self._executor.shutdown(wait=False, cancel_futures=True)
            self._executor = None
        self._pending = None
        self._pending_index = None

    def _load(self, index: int):
        chunk_index = bisect.bisect_right(self._starts, index)-1
        if chunk_index < 0 or index >= self._starts[chunk_index]+self.chunks[chunk_index]["count"]:
            raise IndexError(index)
        if chunk_index != self._chunk_index:
            if chunk_index in self._cache:
                data = self._cache[chunk_index]
                self._cache.move_to_end(chunk_index)
            elif chunk_index == self._pending_index:
                data = self._pending.result()
                self._pending = None
                self._pending_index = None
                self._remember_chunk(chunk_index, data)
            else:
                data = self._read_chunk(chunk_index)
                self._remember_chunk(chunk_index, data)
            self._data = data
            self._chunk_index = chunk_index
        self._prefetch_next(chunk_index)
        return index-self._starts[chunk_index]

    def frame(self, index: int) -> Snapshot:
        if not 0 <= index < len(self):
            raise IndexError(index)
        row = self._load(index)
        data = self._data
        phase = float(data["shaft_phase"][row])
        membership = data["membership"][row].copy() if "membership" in data else None
        return Snapshot(
            float(data["time"][row]), data["position"][row].copy(),
            data["velocity"][row].copy(), data["angle"][row].copy(),
            data["omega"][row].copy(), float(data["energy"][row]),
            float(data["temperature_trans"][row]),
            float(data["temperature_rot"][row]),
            str(data["branch"][row]) or None,
            tuple(int(x) for x in data["occupancy"][row]),
            int(data["event_count"][row]),
            float(data["energy_residual"][row]),
            tuple(ApparatusComponent(**{**x,
                  "points": tuple(tuple(p) for p in x["points"]),
                  "velocity": tuple(x["velocity"])}) for x in
                  json.loads(str(data["apparatus"][row]))),
            None if math.isnan(phase) else phase, membership=membership)

    def seek(self, time: float) -> Snapshot:
        """Return the latest recorded frame at or before a physical time."""
        return self.frame(self.seek_index(time))

    def seek_index(self, time: float) -> int:
        """Index of the latest recorded frame at or before a physical time."""
        if not math.isfinite(time):
            raise ValueError("seek time must be finite")
        if time <= 0:
            return 0
        chunk_index = bisect.bisect_left(
            [c["last_time"] for c in self.chunks], time)
        chunk_index = min(chunk_index, len(self.chunks)-1)
        self._load(self.chunks[chunk_index]["first_frame"])
        times = self._data["time"]
        offset = max(0, bisect.bisect_right(times, time)-1)
        if offset == 0 and time < times[0] and chunk_index > 0:
            previous = self.chunks[chunk_index-1]
            return previous["first_frame"]+previous["count"]-1
        return self.chunks[chunk_index]["first_frame"]+offset

    def ledger(self, index: int) -> dict[str, float]:
        if not 0 <= index < len(self):
            raise IndexError(index)
        row = self._load(index)
        return dict(zip(self.manifest["ledger_columns"],
                        (float(x) for x in self._data["ledger"][row])))

    @property
    def has_instruments(self) -> bool:
        return "instruments" in self.manifest

    def instrument_series(self) -> dict[str, np.ndarray]:
        """Whole-run scalar columns; reads only small columns of each chunk."""
        if self._series is None:
            names = ["time", "temperature_trans", "temperature_rot", "branch",
                     "energy_residual", "event_count", "ledger"]
            if self.has_instruments:
                names += ["i_"+name for name in INSTRUMENT_COLUMNS]
            parts = {name: [] for name in names}
            for chunk in self.chunks:
                with np.load(self.path/chunk["file"], allow_pickle=False) as archive:
                    for name in names:
                        parts[name].append(archive[name])
            self._series = {name.removeprefix("i_"): np.concatenate(values)
                            for name, values in parts.items()}
        return self._series

    @property
    def has_osmosis(self) -> bool:
        return "osmosis" in self.manifest

    def osmosis_frame(self, index: int):
        """The osmosis instruments' frame as the live window showed it (v3)."""
        if not self.has_osmosis:
            raise ValueError("this archive holds no osmosis instruments")
        if not 0 <= index < len(self):
            raise IndexError(index)
        if self._osmosis is None:
            parts = {}
            for chunk in self.chunks:
                with np.load(self.path/chunk["file"], allow_pickle=False) as archive:
                    for name in archive.files:
                        if name.startswith("o_"):
                            parts.setdefault(name[2:], []).append(archive[name])
            columns = {name: np.concatenate(values) for name, values in parts.items()}
            self._osmosis = RecordedOsmosis(columns, self.manifest["osmosis"],
                                            self.manifest["world_metadata"])
        return self._osmosis.frame(index)

    def cycles(self):
        """Completed-cycle summaries and pressure-area comparisons (v2)."""
        if self._cycles is None:
            info = self.manifest.get("instruments", {})
            markers = [CycleMarker(**m) for m in info.get("cycle_markers", [])]
            summaries = [summarize_cycle(a, b) for a, b in zip(markers[:-1], markers[1:])]
            comparisons = [None if c is None else PressureAreaComparison(**c)
                           for c in info.get("cycle_pressure_area", [])]
            self._cycles = (tuple(summaries), tuple(comparisons), tuple(markers))
        return self._cycles

    def instrument_sample(self, index: int) -> InstrumentSample:
        """The live-window instrument reading recorded at a saved frame (v2)."""
        if not self.has_instruments:
            raise ValueError("this archive predates recorded instruments")
        series = self.instrument_series()
        info = self.manifest["instruments"]
        summaries, comparisons, _ = self.cycles()
        completed = int(series["completed_cycles"][index])
        if completed not in self._efficiency:
            self._efficiency[completed] = assess_efficiency(
                summaries[:completed], info["transient_cycles"],
                info["efficiency_min_cycles"])
        ledger = dict(zip(self.manifest["ledger_columns"], series["ledger"][index]))
        metadata = self.manifest["world_metadata"]

        def number(name):
            value = float(series[name][index])
            return None if math.isnan(value) else value
        return InstrumentSample(
            time=float(series["time"][index]), area=number("area"),
            pressure=number("pressure"), pressure_window=info["pressure_window"],
            temperature_trans=float(series["temperature_trans"][index]),
            temperature_rot=float(series["temperature_rot"][index]),
            reservoir_hot=metadata.get("T_hot"), reservoir_cold=metadata.get("T_cold"),
            gas_energy=float(series["gas_energy"][index]),
            piston_energy=float(series["piston_energy"][index]),
            flywheel_energy=float(series["flywheel_energy"][index]),
            spring_energy=float(series["spring_energy"][index]),
            heat_hot=float(ledger["heat_hot"]), heat_cold=float(ledger["heat_cold"]),
            motor_work=float(series["motor_work"][index]),
            load_output=float(ledger["load_output"]),
            first_law_residual=float(series["energy_residual"][index]),
            branch=str(series["branch"][index]) or None,
            event_count=int(series["event_count"][index]),
            max_penetration=float(series["max_penetration"][index]),
            ccd_refinements=int(series["ccd_refinements"][index]),
            ccd_failures=int(series["ccd_failures"][index]),
            cluster_count=int(series["cluster_count"][index]),
            max_cluster_residual=float(series["max_cluster_residual"][index]),
            entropy=number("entropy"), completed_cycles=completed,
            latest_cycle=summaries[completed-1] if completed else None,
            latest_pressure_area=(comparisons[completed-1]
                                  if 0 < completed <= len(comparisons) else None),
            efficiency_report=self._efficiency[completed])

    def _times(self):
        if self._event_times is None:
            with np.load(self.path/self.manifest["event_times_file"],
                         allow_pickle=False) as archive:
                self._event_times = archive["time"]
                self._event_kinds = archive["kind"]
                self._event_participants = archive["participants"]
                if "impulse" in archive.files:
                    self._event_impulse = archive["impulse"]
                if "point" in archive.files:
                    self._event_point = archive["point"]
            participants = [json.loads(str(value)) for value in
                            self._event_participants]
            self._event_a = np.array([ids[0] if ids else -1 for ids in
                                      participants], dtype=np.int32)
            self._event_b = np.array([ids[1] if len(ids) > 1 else -1 for ids in
                                      participants], dtype=np.int32)
            if self._event_impulse is None:
                # Older version-1 archives retain impulses in events.jsonl.
                with (self.path/"events.jsonl").open() as file:
                    impulses = [json.loads(line)["impulse"] for line in file]
                self._event_impulse = np.asarray(impulses,
                                                 dtype=np.float64).reshape(-1, 2)
            if len(self._event_impulse) != len(self._event_times):
                raise ValueError("replay event index and event log disagree")
            if self._event_point is not None and len(self._event_point) != len(self._event_times):
                raise ValueError("replay event points and event index disagree")
        return self._event_times

    def _body_properties(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self._mass is None:
            geometry = self.manifest["geometry"]
            if all(key in geometry for key in ("mass", "inertia", "dynamic")):
                self._mass = np.asarray(geometry["mass"], dtype=np.float64)
                self._inertia = np.asarray(geometry["inertia"], dtype=np.float64)
                self._dynamic = np.asarray(geometry["dynamic"], dtype=bool)
            else:
                # Older version-1 archives retain the deterministic config.
                bodies = load_preset(RunConfig(**self.manifest["config"])).world.bodies
                self._mass = bodies.mass.copy()
                self._inertia = bodies.inertia.copy()
                self._dynamic = bodies.dynamic.copy()
        return self._mass, self._inertia, self._dynamic

    def _can_reconstruct_motion(self) -> bool:
        if self._exact_motion is None:
            self._times()
            if self._exact_disc:
                self._exact_motion = True
            elif self._event_point is None:
                self._exact_motion = False
            else:
                _, inertia, dynamic = self._body_properties()
                active = np.any(self._event_impulse != 0, axis=1)
                angular = np.zeros(len(active), dtype=bool)
                for bodies in (self._event_a, self._event_b):
                    present = bodies >= 0
                    angular[present] |= (np.isfinite(inertia[bodies[present]]) &
                                         dynamic[bodies[present]])
                self._exact_motion = bool(np.all(np.isfinite(
                    self._event_point[active & angular])))
        return self._exact_motion

    def _event_state_at(self, left: Snapshot, end_event: int,
                        time: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Advance particle motion and rotation between recorded impulses."""
        self._times()
        mass, inertia, dynamic = self._body_properties()
        position = left.position.copy()
        velocity = left.velocity.copy()
        angle = left.angle.copy()
        omega = left.omega.copy()
        updated_at = np.full(len(velocity), left.time)
        stop = bisect.bisect_right(self._event_times, time,
                                  left.event_count, end_event)
        for event_index in range(left.event_count, stop):
            event_time = self._event_times[event_index]
            impulse = self._event_impulse[event_index]
            a, b = self._event_a[event_index], self._event_b[event_index]
            if a >= 0 and dynamic[a]:
                dt = event_time - updated_at[a]
                position[a] += velocity[a] * dt
                angle[a] = (angle[a] + omega[a] * dt + math.pi) % (2 * math.pi) - math.pi
                updated_at[a] = event_time
                jv = -impulse
                velocity[a] += jv / mass[a]
                if math.isfinite(inertia[a]):
                    point = self._event_point[event_index]
                    arm = point - position[a]
                    omega[a] += (arm[0] * jv[1] - arm[1] * jv[0]) / inertia[a]
            if b >= 0 and dynamic[b]:
                dt = event_time - updated_at[b]
                position[b] += velocity[b] * dt
                angle[b] = (angle[b] + omega[b] * dt + math.pi) % (2 * math.pi) - math.pi
                updated_at[b] = event_time
                jv = impulse
                velocity[b] += jv / mass[b]
                if math.isfinite(inertia[b]):
                    point = self._event_point[event_index]
                    arm = point - position[b]
                    omega[b] += (arm[0] * jv[1] - arm[1] * jv[0]) / inertia[b]
        position[dynamic] += (velocity[dynamic] *
                              (time - updated_at[dynamic])[:, None])
        angle[dynamic] = np.remainder(
            angle[dynamic] + omega[dynamic] * (time - updated_at[dynamic]) + math.pi,
            2 * math.pi) - math.pi
        return position, velocity, angle, omega

    def sample(self, time: float) -> tuple[Snapshot, int]:
        """Reconstruct event-aware motion when contact points are available.

        Scalar statistics remain at the preceding saved timestamp. Earlier
        archives without angular contact points use visual interpolation for
        non-disc shapes.
        """
        index = self.seek_index(time)
        left = self.frame(index)
        if index+1 >= len(self) or time <= left.time:
            return left, index
        right = self.frame(index+1)
        alpha = (time-left.time)/(right.time-left.time)
        exact_motion = self._can_reconstruct_motion()
        if exact_motion:
            position, velocity, angle, omega = self._event_state_at(
                left, right.event_count, time)
        else:
            position = left.position+(right.position-left.position)*alpha
            angle = left.angle+(right.angle-left.angle)*alpha
        phase = (None if left.shaft_phase is None or right.shaft_phase is None
                 else left.shaft_phase+(right.shaft_phase-left.shaft_phase)*alpha)
        branch = left.branch
        if right.branch != left.branch and right.event_count > left.event_count:
            self._times()
            start, end = left.event_count, right.event_count
            transitions = self._event_times[start:end][
                self._event_kinds[start:end] == "branch_transition"]
            if len(transitions) and time >= transitions[0]:
                branch = right.branch
        branch_switched = branch != left.branch
        membership = left.membership
        if membership is not None and right.event_count > left.event_count:
            # Discs enter and leave wells at logged step events.
            self._times()
            start = left.event_count
            stop = bisect.bisect_right(self._event_times, time, start, right.event_count)
            steps = np.flatnonzero(self._event_kinds[start:stop] == "step")
            if len(steps):
                membership = membership.copy()
                discs = self._event_a[start + steps]
                membership[discs] = right.membership[discs]
        apparatus = _interpolate_apparatus(left, right, alpha, branch_switched)
        # Earlier archives lack angular contact points. Keep their visual
        # approximation until they are regenerated with this information.
        if not exact_motion and right.event_count > left.event_count:
            self._times()
            velocity = left.velocity.copy()
            omega = left.omega.copy()
            last_event = np.full(len(velocity), math.inf)
            for event_index in range(left.event_count,right.event_count):
                for body in (self._event_a[event_index],self._event_b[event_index]):
                    if body < 0:
                        continue
                    last_event[body] = self._event_times[event_index]
            settled = last_event <= time
            velocity[settled] = right.velocity[settled]
            omega[settled] = right.omega[settled]
        elif not exact_motion:
            velocity = left.velocity
            omega = left.omega
        return replace(left,time=time,position=position,velocity=velocity,
                       angle=angle,omega=omega,shaft_phase=phase,
                       branch=branch,apparatus=apparatus,membership=membership), index
