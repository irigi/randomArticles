"""Cycle work from the causal piston-impulse pressure instrument."""

from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass

from .cycles import CycleMarker, summarize_cycle


@dataclass(frozen=True)
class PressureAreaComparison:
    """Work by gas; the bound covers window smoothing and cycle-end clipping."""

    integrated_work_by_gas: float
    event_work_by_gas: float
    smoothing_error_bound: float
    window: float

    @property
    def difference(self) -> float:
        return self.integrated_work_by_gas-self.event_work_by_gas


def piston_wall_index(simulation) -> int | None:
    return next((index for index, wall in enumerate(simulation.world.walls)
                 if wall.name == "cam_piston"), None)


def piston_impulses(simulation, events, first_index: int):
    """Yield (global event index, time, |normal impulse|) for piston contacts."""
    piston = piston_wall_index(simulation)
    for offset, event in enumerate(events):
        if (event.metadata or {}).get("boundary") == piston:
            yield first_index+offset, event.time, abs(event.impulse[0])


def compare_pressure_area(simulation, start: CycleMarker, end: CycleMarker,
                          window: float = .25) -> PressureAreaComparison | None:
    """Integrate the trailing pressure window over the retained event history."""
    if window <= 0:
        raise ValueError("pressure window must be positive")
    first = bisect_left(simulation.events, start.time-window,
                        hi=end.event_count, key=lambda event: event.time)
    impulses = piston_impulses(simulation, simulation.events[first:end.event_count], first)
    return pressure_area_from_impulses(simulation, start, end, impulses, window)


def pressure_area_from_impulses(simulation, start: CycleMarker, end: CycleMarker,
                                impulses, window: float = .25
                                ) -> PressureAreaComparison | None:
    """Integrate the fixed-width trailing pressure window exactly for a controlled cam.

    ``impulses`` holds (global event index, time, |impulse|) piston contacts
    covering at least [start.time-window, end.time]. Each impulse contributes
    constant pressure for ``window`` seconds. Its integral against dA is the
    exact cam-area difference over the portion of that interval inside this
    cycle. The error bound sums absolute per-impulse changes caused by
    smoothing and clipping; it is deterministic, not a statistical confidence
    interval. Free-shaft histories need stored phase trajectories before this
    estimator can be used there.
    """
    if window <= 0:
        raise ValueError("pressure window must be positive")
    mechanism = simulation.world.mechanism
    shaft = getattr(mechanism, "shaft", None)
    cam = getattr(mechanism, "cam", None)
    if shaft is None or cam is None or shaft.prescribed_omega is None:
        return None
    if piston_wall_index(simulation) is None:
        return None
    omega = shaft.prescribed_omega

    def phase(time: float) -> float:
        return start.phase+omega*(time-start.time)

    integrated = 0.
    bound = 0.
    for index, time, magnitude in impulses:
        if index >= end.event_count or time < start.time-window:
            continue
        impulse = magnitude/cam.height
        left = max(start.time, time)
        right = min(end.time, time+window)
        smooth = (impulse*(cam.area(phase(right))-cam.area(phase(left)))/window
                  if right > left else 0.)
        instantaneous = (impulse*cam.area_derivative(phase(time))*omega
                         if index >= start.event_count else 0.)
        integrated += smooth
        bound += abs(smooth-instantaneous)
    event_work = -summarize_cycle(start, end).piston_work_on_gas
    return PressureAreaComparison(integrated, event_work, bound, window)
