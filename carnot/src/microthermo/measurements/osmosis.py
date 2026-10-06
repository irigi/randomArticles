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

    def __init__(self, history_limit: int = 3000, history_spacing: float = .25,
                 force_window: float = 10., average_after: float = 20.):
        self.history: deque[OsmosisSample] = deque(maxlen=history_limit)
        self.history_spacing = history_spacing
        self.force_window = force_window
        self.average_after = average_after
        self._impulses: deque[tuple[float, float]] = deque()
        self._height = 1.
        self._reset_averages()

    def observe(self, simulation, snapshot=None) -> OsmosisFrame:
        snapshot = simulation.snapshot() if snapshot is None else snapshot
        meta = simulation.world.metadata
        t = snapshot.time
        if self.history and t < self.history[-1].time:      # restored or reset
            self.history.clear()
            self._impulses.clear()
            self._reset_averages()
        hosts = int(meta["hosts"])
        temperature = float(meta["temperature"])
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
        if t >= self.average_after:
            if self._impulse_start is None:
                self._impulse_start = (t, disc_impulse, host_impulse)
            self._impulse_now = (t, disc_impulse, host_impulse)
        if not self.history or t - self.history[-1].time >= self.history_spacing:
            self.history.append(sample)
            if t >= self.average_after:
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
        """Means over the samples after ``average_after``. The osmotic pressure
        is the net post impulse over that whole interval per membrane height,
        with its split by contacting body."""
        if not self._count:
            return {}
        left, right, bound, free_left = self._sums/self._count
        averages = {"samples": self._count, "left": left, "right": right, "bound": bound,
                    "free_left": free_left}
        (t0, d0, h0), (t1, d1, h1) = self._impulse_start, self._impulse_now
        if t1 > t0:
            span = (t1 - t0)*self._height
            averages["osmotic_pressure"] = (h1 - h0 + d1 - d0)/span
            averages["host_contact_part"] = (h1 - h0)/span
            averages["disc_contact_part"] = (d1 - d0)/span
            averages["duration"] = t1 - t0
        return averages
