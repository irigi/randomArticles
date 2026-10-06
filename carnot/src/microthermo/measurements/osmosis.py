"""Passive osmosis instruments: counts, chemical potentials and membrane force.

Every quantity comes from the bodies' positions, the kernel's well
membership and the impulse delivered to the membrane posts. The ideal-point
reference (experiments.osmosis.OsmosisDesign.reference) is shown beside the
measurements and is never fed back into the dynamics.

The osmotic pressure is the total net x-force on the membrane posts per
membrane height. It is also split by the body in each post contact. That
split is a diagnostic, not a pressure of its own: discs push hosts against
the membrane (fewer discs fit between a host and the posts), so with more
discs the hosts' part grows and the discs' part turns negative while the
total stays put.

Averages (side counts, bound discs, the occupancy histogram and the
pressure) start ``average_after`` into the run and start again that long
after every change of the wall temperature, so they always describe one
temperature.

Chemical potentials are reported as mu/T = ln(c) relative to unit density:
c_R on the right, c_free on the left, and n_bound/(a e^{eps/T}) per well.
At equilibrium of ideal points all three agree.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math

import numpy as np

from ..experiments.osmosis import ideal_reference


@dataclass(frozen=True)
class OsmosisSample:
    time: float
    left: int
    right: int
    bound: int
    free_left: int
    per_host: tuple[int, ...]
    mu_right: float
    mu_free_left: float
    mu_bound: float
    membrane_force: float | None       # net x-force on the posts over the window
    osmotic_pressure: float | None     # that force per membrane height
    host_contact_part: float | None    # from host contacts, per height
    disc_contact_part: float | None    # from disc contacts (left minus right), per height
    temperature: float


@dataclass(frozen=True)
class OsmosisFrame:
    current: OsmosisSample
    history: tuple[OsmosisSample, ...]
    reference: dict
    occupancy_histogram: np.ndarray    # time-sampled counts of n per host, n = 0..
    averages: dict


def _log(value: float) -> float:
    return math.log(value) if value > 0 else math.nan


class OsmosisInstruments:
    """Observe an osmosis world; ``observe`` never changes the simulation."""

    def __init__(self, history_limit: int = 4000, history_spacing: float = .25,
                 force_window: float = 10., average_after: float = 20.):
        self.history: deque[OsmosisSample] = deque(maxlen=history_limit)
        self.history_spacing = history_spacing
        self.force_window = force_window
        self.average_after = average_after
        self._impulses: deque[tuple[float, float]] = deque()
        self._height = 1.
        self._temperature = None
        self.segment = 0                   # averaging restarts so far
        self._average_from = average_after
        self._reset_averages()

    def observe(self, simulation, snapshot=None) -> OsmosisFrame:
        snapshot = simulation.snapshot() if snapshot is None else snapshot
        meta = simulation.world.metadata
        t = snapshot.time
        if self.history and t < self.history[-1].time:      # restored or reset
            self.history.clear()
            self._impulses.clear()
            self._reset_averages()
            self._temperature = None
            self._average_from = self.average_after
        hosts = int(meta["hosts"])
        temperature = float(meta["temperature"])
        if self._temperature is not None and temperature != self._temperature:
            self._reset_averages()
            self.segment += 1
            self._average_from = t + self.average_after
        self._temperature = temperature
        eps = float(meta["ring"]["well_depth"])
        x = snapshot.position[np.asarray(simulation.bh) < 0, 0]
        left = int(np.sum(x < float(meta["membrane_x"])))
        right = len(x) - left
        per_host = tuple(int(v) for v in simulation.occupancy())
        bound = sum(per_host)
        free_left = left - bound
        weight = float(meta["well_area"])*math.exp(eps/temperature)
        mu_right = _log(right/float(meta["accessible_area_right"]))
        mu_free = _log(free_left/float(meta["free_left_area"]))
        mu_bound = _log(bound/(hosts*weight))
        # Net x-impulse on the posts (positive when the left side pushes
        # harder), split by the moving body of each contact.
        disc_impulse = float(simulation.post_impulses("discs")[:, 0].sum())
        host_impulse = float(simulation.post_impulses("hosts")[:, 0].sum())
        self._height = simulation.world.walls[0].length
        self._impulses.append((t, disc_impulse, host_impulse))
        while len(self._impulses) > 2 and self._impulses[1][0] <= t - self.force_window:
            self._impulses.popleft()
        t0, d0, h0 = self._impulses[0]
        if t > t0:
            disc_rate, host_rate = (disc_impulse - d0)/(t - t0), (host_impulse - h0)/(t - t0)
            force = disc_rate + host_rate
            pressure, host_part = force/self._height, host_rate/self._height
            disc_part = disc_rate/self._height
        else:
            force = pressure = host_part = disc_part = None
        sample = OsmosisSample(t, left, right, bound, free_left, per_host, mu_right, mu_free,
                               mu_bound, force, pressure, host_part, disc_part, temperature)
        if t >= self._average_from:
            if self._impulse_start is None:
                self._impulse_start = (t, disc_impulse, host_impulse)
            self._impulse_now = (t, disc_impulse, host_impulse)
        if not self.history or t - self.history[-1].time >= self.history_spacing:
            self.history.append(sample)
            if t >= self._average_from:
                self._accumulate(sample)
        reference = ideal_reference(
            float(meta["accessible_area_right"]), float(meta["free_left_area"]),
            float(meta["well_area"]), eps, hosts, float(meta["host_accessible_area"]),
            len(x), temperature)
        return OsmosisFrame(sample, tuple(self.history), reference,
                            self._histogram.copy(), self.averages())

    def _reset_averages(self) -> None:
        self._histogram = np.zeros(1, np.int64)
        self._sums = None
        self._count = 0
        self._impulse_start = self._impulse_now = None

    def _accumulate(self, sample: OsmosisSample) -> None:
        top = max(sample.per_host, default=0)
        if top >= len(self._histogram):
            self._histogram = np.pad(self._histogram, (0, top + 1 - len(self._histogram)))
        for n in sample.per_host:
            self._histogram[n] += 1
        values = np.array([sample.left, sample.right, sample.bound, sample.free_left])
        self._sums = values if self._sums is None else self._sums + values
        self._count += 1

    def averages(self) -> dict:
        """Means over the samples since averaging last (re)started. The osmotic
        pressure is the net post impulse over that interval per membrane
        height, with its split by contacting body."""
        if not self._count:
            return {}
        left, right, bound, free_left = self._sums/self._count
        averages = {"samples": self._count, "left": left, "right": right, "bound": bound,
                    "free_left": free_left, "since": self._impulse_start[0]}
        averages.update(_pressure_averages(self._impulse_start, self._impulse_now,
                                           self._height))
        return averages


def _pressure_averages(start, now, height: float) -> dict:
    (t0, d0, h0), (t1, d1, h1) = start, now
    if not t1 > t0:
        return {}
    span = (t1 - t0)*height
    return {"osmotic_pressure": (h1 - h0 + d1 - d0)/span,
            "host_contact_part": (h1 - h0)/span, "disc_contact_part": (d1 - d0)/span,
            "duration": t1 - t0}


# Per-frame replay columns (archive version 3), stored with the prefix "o_".
# ``in_history`` marks the samples the live instruments kept for plotting
# and averaging; with the raw post impulses it lets a replay rebuild the
# exact OsmosisFrame of any recorded frame.
SAMPLE_COLUMNS = {
    "time": np.float64, "left": np.int64, "right": np.int64, "bound": np.int64,
    "free_left": np.int64, "mu_right": np.float64, "mu_free_left": np.float64,
    "mu_bound": np.float64, "membrane_force": np.float64,
    "osmotic_pressure": np.float64, "host_contact_part": np.float64,
    "disc_contact_part": np.float64, "temperature": np.float64,
}
RAW_COLUMNS = {"disc_impulse": np.float64, "host_impulse": np.float64,
               "in_history": np.bool_, "segment": np.int64, "average_from": np.float64}


class OsmosisRecorder:
    """Run OsmosisInstruments while precalculating and collect replay columns."""

    def __init__(self, instruments: OsmosisInstruments | None = None):
        self.instruments = OsmosisInstruments() if instruments is None else instruments
        self.rows: list[tuple] = []

    def settings(self, simulation) -> dict:
        i = self.instruments
        return {"history_limit": i.history.maxlen, "history_spacing": i.history_spacing,
                "force_window": i.force_window, "average_after": i.average_after,
                "height": simulation.world.walls[0].length,
                "discs": int(np.sum(np.asarray(simulation.bh) < 0))}

    def observe(self, simulation, snapshot) -> OsmosisFrame:
        before = self.instruments.history[-1] if self.instruments.history else None
        frame = self.instruments.observe(simulation, snapshot)
        kept = bool(self.instruments.history) and self.instruments.history[-1] is not before
        t, disc, host = self.instruments._impulses[-1]
        self.rows.append((frame.current, disc, host, kept, self.instruments.segment,
                          self.instruments._average_from))
        return frame

    def take_columns(self) -> dict[str, np.ndarray]:
        columns = {}
        for name, dtype in SAMPLE_COLUMNS.items():
            values = [getattr(row[0], name) for row in self.rows]
            columns["o_"+name] = np.array([math.nan if v is None else v for v in values],
                                          dtype=dtype)
        for k, (name, dtype) in enumerate(RAW_COLUMNS.items(), start=1):
            columns["o_"+name] = np.array([row[k] for row in self.rows], dtype=dtype)
        columns["o_per_host"] = np.array([row[0].per_host for row in self.rows],
                                         dtype=np.int64).reshape(len(self.rows), -1)
        self.rows.clear()
        return columns


class RecordedOsmosis:
    """Rebuild the live OsmosisFrame at any frame of a version-3 replay."""

    def __init__(self, columns: dict[str, np.ndarray], settings: dict, metadata: dict):
        self.columns = columns
        self.settings = settings
        self.metadata = metadata
        self.kept = np.flatnonzero(columns["in_history"])
        self._samples: dict[int, OsmosisSample] = {}
        time, segment = columns["time"], columns["segment"]
        averaged = np.flatnonzero(columns["in_history"] & (time >= columns["average_from"]))
        self.averaged = averaged
        self._averaged_segment = segment[averaged]
        values = np.stack([columns[n][averaged] for n in ("left", "right", "bound",
                                                          "free_left")], axis=1)
        self._sums = np.cumsum(values, axis=0)
        per_host = columns["per_host"][averaged]
        top = int(per_host.max(initial=0)) + 1
        counts = np.zeros((len(averaged), top), np.int64)
        for n in range(top):
            counts[:, n] = np.sum(per_host == n, axis=1)
        self._histograms = np.cumsum(counts, axis=0)
        # First frame of each averaging segment at which averaging is on.
        on = np.flatnonzero(time >= columns["average_from"])
        self._impulse_start = {int(segment[k]): int(k) for k in on[::-1]}

    def sample(self, index: int) -> OsmosisSample:
        if index not in self._samples:
            c = self.columns

            def number(name):
                value = float(c[name][index])
                return None if math.isnan(value) else value
            self._samples[index] = OsmosisSample(
                float(c["time"][index]), int(c["left"][index]), int(c["right"][index]),
                int(c["bound"][index]), int(c["free_left"][index]),
                tuple(int(v) for v in c["per_host"][index]), float(c["mu_right"][index]),
                float(c["mu_free_left"][index]), float(c["mu_bound"][index]),
                number("membrane_force"), number("osmotic_pressure"),
                number("host_contact_part"), number("disc_contact_part"),
                float(c["temperature"][index]))
        return self._samples[index]

    def frame(self, index: int) -> OsmosisFrame:
        c, settings, meta = self.columns, self.settings, self.metadata
        stop = int(np.searchsorted(self.kept, index, side="right"))
        history = tuple(self.sample(int(k)) for k in
                        self.kept[max(0, stop - settings["history_limit"]):stop])
        end = int(np.searchsorted(self.averaged, index, side="right"))
        segment = int(c["segment"][index])
        begin = int(np.searchsorted(self._averaged_segment[:end], segment, side="left"))
        count = end - begin
        averages = {}
        histogram = np.zeros(1, np.int64)
        if count:
            def window(cumulative):
                return cumulative[end - 1] - (cumulative[begin - 1] if begin else 0)
            left, right, bound, free_left = window(self._sums)/count
            k0 = self._impulse_start[segment]
            averages = {"samples": count, "left": left, "right": right, "bound": bound,
                        "free_left": free_left, "since": float(c["time"][k0])}
            histogram = window(self._histograms)
            histogram = histogram[:max(1, int(np.flatnonzero(histogram).max(initial=0)) + 1)]

            def point(k):
                return (float(c["time"][k]), float(c["disc_impulse"][k]),
                        float(c["host_impulse"][k]))
            averages.update(_pressure_averages(point(k0), point(index), settings["height"]))
        temperature = float(c["temperature"][index])
        reference = ideal_reference(
            float(meta["accessible_area_right"]), float(meta["free_left_area"]),
            float(meta["well_area"]), float(meta["ring"]["well_depth"]), int(meta["hosts"]),
            float(meta["host_accessible_area"]), int(settings["discs"]), temperature)
        return OsmosisFrame(self.sample(index), history, reference, histogram, averages)
