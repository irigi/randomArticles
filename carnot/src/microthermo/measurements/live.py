"""Passive, bounded live measurements for the desktop laboratory."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from .cycles import CycleSummary, EfficiencyReport, assess_efficiency, summarize_cycle
from .pressure_area import PressureAreaComparison, compare_pressure_area
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
    completed_cycles: int
    latest_cycle: CycleSummary | None
    latest_pressure_area: PressureAreaComparison | None
    efficiency_report: EfficiencyReport


@dataclass(frozen=True)
class InstrumentFrame:
    current: InstrumentSample
    history: tuple[InstrumentSample, ...]
    ideal_reference: IdealCarnotReference | None


class LiveInstruments:
    """Read-only subscriber; observations never advance or modify physics."""

    def __init__(self, history_limit: int = 600, pressure_window: float = .25,
                 transient_cycles: int = 2, efficiency_min_cycles: int = 8):
        if history_limit < 2 or pressure_window <= 0:
            raise ValueError("require at least two samples and positive pressure window")
        if transient_cycles < 0 or efficiency_min_cycles < 4:
            raise ValueError("invalid efficiency cycle gate")
        self.history = deque(maxlen=history_limit)
        self.impulses = deque()
        self.pressure_window = pressure_window
        self.event_cursor = 0
        self._comparison_markers = None
        self._latest_pressure_area = None
        self._ideal_reference = None
        self._reference_ready = False
        self.transient_cycles=transient_cycles
        self.efficiency_min_cycles=efficiency_min_cycles
        self._efficiency_markers=None
        self._efficiency_report=assess_efficiency((),transient_cycles,
                                                  efficiency_min_cycles)

    def observe(self, simulation, snapshot=None) -> InstrumentFrame:
        snapshot = simulation.snapshot() if snapshot is None else snapshot
        if not self._reference_ready:
            self._ideal_reference=ideal_carnot_reference(simulation)
            self._reference_ready=True
        if self.event_cursor > len(simulation.events):
            self.event_cursor = 0
            self.impulses.clear()
            self.history.clear()
            self._comparison_markers = None
            self._latest_pressure_area = None
            self._efficiency_markers=None
        mechanism = simulation.world.mechanism
        cam = getattr(mechanism,"cam",None)
        piston_wall = next((index for index,wall in
                            enumerate(simulation.world.walls)
                            if wall.name == "cam_piston"),None)
        if cam is not None and piston_wall is not None:
            for event in simulation.events[self.event_cursor:]:
                if (event.metadata or {}).get("boundary") == piston_wall:
                    self.impulses.append((event.time,abs(event.impulse[0])))
        self.event_cursor = len(simulation.events)
        time = snapshot.time
        while self.impulses and self.impulses[0][0] <= time-self.pressure_window:
            self.impulses.popleft()
        # Zero-pad the pre-run part of the fixed trailing window. This is the
        # same pressure definition integrated by compare_pressure_area().
        window = self.pressure_window
        pressure = (sum(impulse for _,impulse in self.impulses)/
                    (cam.height*window) if cam is not None and time > 0 else None)
        area = cam.area(mechanism.shaft.phi) if cam is not None else None
        parts = {part.name:part for part in snapshot.apparatus}
        ledger = simulation.ledger
        markers=simulation.cycle_markers
        latest_cycle=(summarize_cycle(markers[-2],markers[-1])
                      if len(markers) >= 2 else None)
        comparison_markers=tuple(markers[-2:]) if len(markers) >= 2 else None
        if comparison_markers != self._comparison_markers:
            self._comparison_markers=comparison_markers
            self._latest_pressure_area=(compare_pressure_area(
                simulation,*comparison_markers,self.pressure_window)
                if comparison_markers else None)
        latest_pressure_area=self._latest_pressure_area
        marker_key=tuple(markers)
        if marker_key != self._efficiency_markers:
            self._efficiency_markers=marker_key
            summaries=[summarize_cycle(start,end) for start,end in
                       zip(markers[:-1],markers[1:])]
            self._efficiency_report=assess_efficiency(
                summaries,self.transient_cycles,self.efficiency_min_cycles)
        sample = InstrumentSample(
            time=time,area=area,pressure=pressure,pressure_window=window,
            temperature_trans=snapshot.translational_temperature,
            temperature_rot=snapshot.rotational_temperature,
            reservoir_hot=simulation.world.metadata.get("T_hot"),
            reservoir_cold=simulation.world.metadata.get("T_cold"),
            gas_energy=simulation.world.bodies.kinetic_energy(),
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
            completed_cycles=max(0,len(markers)-1),latest_cycle=latest_cycle,
            latest_pressure_area=latest_pressure_area,
            efficiency_report=self._efficiency_report)
        self.history.append(sample)
        return InstrumentFrame(sample,tuple(self.history),self._ideal_reference)
