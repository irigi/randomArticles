"""Passive, bounded live measurements for the desktop laboratory."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math

import numpy as np
from .cycles import CycleMarker, CycleSummary, EfficiencyReport, assess_efficiency, summarize_cycle
from .pressure_area import (PressureAreaComparison, piston_impulses,
                            pressure_area_from_impulses)
from .ideal import IdealCarnotReference, ideal_carnot_reference


@dataclass(frozen=True)
class InstrumentSample:
    time: float
    area: float | None
    pressure: float | None
    pressure_window: float
    temperature_trans: float
    temperature_rot: float
    reservoir_hot: float | None
    reservoir_cold: float | None
    gas_energy: float
    piston_energy: float
    flywheel_energy: float
    spring_energy: float
    heat_hot: float
    heat_cold: float
    motor_work: float
    load_output: float
    first_law_residual: float
    branch: str | None
    event_count: int
    max_penetration: float
    ccd_refinements: int
    ccd_failures: int
    cluster_count: int
    max_cluster_residual: float
    # Ideal-gas entropy per particle (k_B = 1) from area and gas energy:
    # ln A + (f/2) ln T with T = 2E/(fN). An estimate, not a measurement.
    entropy: float | None
    completed_cycles: int
    latest_cycle: CycleSummary | None
    latest_pressure_area: PressureAreaComparison | None
    efficiency_report: EfficiencyReport


@dataclass(frozen=True)
class InstrumentFrame:
    current: InstrumentSample
    history: tuple[InstrumentSample, ...]
    ideal_reference: IdealCarnotReference | None
    # Completed cycles in order, and the marker that opened the current one.
    cycles: tuple[CycleSummary, ...] = ()
    cycle_start: CycleMarker | None = None
    particles: int = 0
    degrees_of_freedom: int = 2


def gas_degrees_of_freedom(simulation) -> int:
    """Two for spinless discs, three when every body can rotate."""
    return 3 if np.all(np.isfinite(simulation.world.bodies.inertia)) else 2


def ideal_gas_entropy(area: float | None, gas_energy: float, particles: int,
                      degrees_of_freedom: int) -> float | None:
    if area is None or area <= 0 or gas_energy <= 0 or particles <= 0:
        return None
    temperature = 2*gas_energy/(degrees_of_freedom*particles)
    return math.log(area)+degrees_of_freedom/2*math.log(temperature)


class LiveInstruments:
    """Read-only subscriber; observations never advance or modify physics.

    Events reach the observer either from ``simulation.events`` (the live
    window keeps the full list) or explicitly through ``observe(events=...)``
    when a caller drains the simulation's event list, as precalculation does.
    History samples are kept at least ``history_spacing`` apart in physical
    time so a bounded history can span whole cycles.
    """

    def __init__(self, history_limit: int = 600, pressure_window: float = .25,
                 transient_cycles: int = 2, efficiency_min_cycles: int = 8,
                 history_spacing: float = 0.):
        if history_limit < 2 or pressure_window <= 0 or history_spacing < 0:
            raise ValueError("require at least two samples and positive pressure window")
        if transient_cycles < 0 or efficiency_min_cycles < 4:
            raise ValueError("invalid efficiency cycle gate")
        self.history = deque(maxlen=history_limit)
        self.history_spacing = history_spacing
        # (global event index, time, |impulse|) piston contacts still needed
        # for the trailing pressure window or the open cycle's P dA integral.
        self.impulses = deque()
        # (time, |impulse|) inside the trailing pressure window only.
        self.window_impulses = deque()
        self.pressure_window = pressure_window
        self.event_cursor = 0
        self._comparison_markers = None
        self._latest_pressure_area = None
        self._ideal_reference = None
        self._reference_ready = False
        self.transient_cycles=transient_cycles
        self.efficiency_min_cycles=efficiency_min_cycles
        self._efficiency_markers=None
        self._summaries=()
        self._efficiency_report=assess_efficiency((),transient_cycles,
                                                  efficiency_min_cycles)
        self.cycle_pressure_area: list[PressureAreaComparison | None] = []

    def _reset(self):
        self.event_cursor = 0
        self.impulses.clear()
        self.window_impulses.clear()
        self.history.clear()
        self._comparison_markers = None
        self._latest_pressure_area = None
        self._efficiency_markers=None
        self.cycle_pressure_area.clear()

    def observe(self, simulation, snapshot=None, events=None) -> InstrumentFrame:
        snapshot = simulation.snapshot() if snapshot is None else snapshot
        if not self._reference_ready:
            self._ideal_reference=ideal_carnot_reference(simulation)
            self._dof=gas_degrees_of_freedom(simulation)
            self._reference_ready=True
        if events is None:
            if self.event_cursor > len(simulation.events):
                self._reset()
            events = simulation.events[self.event_cursor:]
            self.event_cursor = len(simulation.events)
        mechanism = simulation.world.mechanism
        cam = getattr(mechanism,"cam",None)
        if cam is not None:
            new = list(piston_impulses(simulation,events,
                                       snapshot.event_count-len(events)))
            self.impulses.extend(new)
            self.window_impulses.extend((t,impulse) for _,t,impulse in new)
        time = snapshot.time
        window = self.pressure_window
        markers=simulation.cycle_markers
        latest_cycle=(summarize_cycle(markers[-2],markers[-1])
                      if len(markers) >= 2 else None)
        comparison_markers=tuple(markers[-2:]) if len(markers) >= 2 else None
        if comparison_markers != self._comparison_markers:
            self._comparison_markers=comparison_markers
            self._latest_pressure_area=(pressure_area_from_impulses(
                simulation,*comparison_markers,tuple(self.impulses),window)
                if comparison_markers else None)
            if comparison_markers:
                self.cycle_pressure_area.append(self._latest_pressure_area)
        latest_pressure_area=self._latest_pressure_area
        # Keep what the open cycle's integral and the trailing window need.
        keep_from = min(time, markers[-1].time if markers else time)-window
        while self.impulses and self.impulses[0][1] < keep_from:
            self.impulses.popleft()
        while self.window_impulses and self.window_impulses[0][0] <= time-window:
            self.window_impulses.popleft()
        # Zero-pad the pre-run part of the fixed trailing window. This is the
        # same pressure definition integrated by compare_pressure_area().
        pressure = (math.fsum(impulse for _,impulse in self.window_impulses)/
                    (cam.height*window) if cam is not None and time > 0 else None)
        area = cam.area(mechanism.shaft.phi) if cam is not None else None
        parts = {part.name:part for part in snapshot.apparatus}
        ledger = simulation.ledger
        marker_key=tuple(markers)
        if marker_key != self._efficiency_markers:
            self._efficiency_markers=marker_key
            summaries=tuple(summarize_cycle(start,end) for start,end in
                            zip(markers[:-1],markers[1:]))
            self._summaries=summaries
            self._efficiency_report=assess_efficiency(
                summaries,self.transient_cycles,self.efficiency_min_cycles)
        gas_energy=simulation.world.bodies.kinetic_energy()
        sample = InstrumentSample(
            time=time,area=area,pressure=pressure,pressure_window=window,
            temperature_trans=snapshot.translational_temperature,
            temperature_rot=snapshot.rotational_temperature,
            reservoir_hot=simulation.world.metadata.get("T_hot"),
            reservoir_cold=simulation.world.metadata.get("T_cold"),
            gas_energy=gas_energy,
            piston_energy=parts.get("piston").kinetic_energy if "piston" in parts else 0.,
            flywheel_energy=parts.get("flywheel").kinetic_energy if "flywheel" in parts else 0.,
            spring_energy=parts.get("shaft_spring").potential_energy if "shaft_spring" in parts else 0.,
            heat_hot=ledger.heat_hot.value,heat_cold=ledger.heat_cold.value,
            motor_work=getattr(mechanism,"motor_work",0.),
            load_output=ledger.load_output.value,
            first_law_residual=snapshot.energy_residual,branch=snapshot.branch,
            event_count=snapshot.event_count,max_penetration=simulation.max_penetration,
            ccd_refinements=simulation.ccd_refinements,
            ccd_failures=simulation.ccd_failures,cluster_count=simulation.cluster_count,
            max_cluster_residual=simulation.max_cluster_residual,
            entropy=ideal_gas_entropy(area,gas_energy,simulation.world.bodies.n,self._dof),
            completed_cycles=max(0,len(markers)-1),latest_cycle=latest_cycle,
            latest_pressure_area=latest_pressure_area,
            efficiency_report=self._efficiency_report)
        if (not self.history or time < self.history[-1].time or
                time-self.history[-1].time >= self.history_spacing):
            self.history.append(sample)
        elif self.history[-1].branch != sample.branch:
            self.history.append(sample)
        return InstrumentFrame(sample,tuple(self.history),self._ideal_reference,
                               self._summaries,markers[-1] if markers else None,
                               simulation.world.bodies.n,self._dof)
