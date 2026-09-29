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


def compare_pressure_area(simulation, start: CycleMarker, end: CycleMarker,
                          window: float = .25) -> PressureAreaComparison | None:
    """Integrate the fixed-width trailing pressure window exactly for a controlled cam.

    Each impulse contributes constant pressure for ``window`` seconds. Its
    integral against dA is the exact cam-area difference over the portion of
    that interval inside this cycle. The error bound sums absolute per-impulse
    changes caused by smoothing and clipping; it is deterministic, not a
    statistical confidence interval. Free-shaft histories need stored phase
    trajectories before this estimator can be used there.
    """
    if window <= 0:
        raise ValueError("pressure window must be positive")
    mechanism = simulation.world.mechanism
    shaft = getattr(mechanism, "shaft", None)
    cam = getattr(mechanism, "cam", None)
    if shaft is None or cam is None or shaft.prescribed_omega is None:
        return None
    piston = next((index for index, wall in enumerate(simulation.world.walls)
                   if wall.name == "cam_piston"), None)
    if piston is None:
        return None
    omega = shaft.prescribed_omega

    def phase(time: float) -> float:
        return start.phase+omega*(time-start.time)

    integrated = 0.
    bound = 0.
    first = bisect_left(simulation.events, start.time-window,
                        hi=end.event_count, key=lambda event: event.time)
    for index in range(first, end.event_count):
        event = simulation.events[index]
        if (event.metadata or {}).get("boundary") != piston:
            continue
        impulse = abs(event.impulse[0])/cam.height
        left = max(start.time, event.time)
        right = min(end.time, event.time+window)
        smooth = (impulse*(cam.area(phase(right))-cam.area(phase(left)))/window
                  if right > left else 0.)
        instantaneous = (impulse*cam.area_derivative(phase(event.time))*omega
                         if index >= start.event_count else 0.)
        integrated += smooth
        bound += abs(smooth-instantaneous)
    event_work = -summarize_cycle(start, end).piston_work_on_gas
    return PressureAreaComparison(integrated, event_work, bound, window)
