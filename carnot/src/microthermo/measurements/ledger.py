from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class CompensatedCounter:
    value: float = 0.0
    correction: float = 0.0

    def add(self, x: float) -> None:
        y = x-self.correction
        t = self.value+y
        self.correction = (t-self.value)-y
        self.value = t


@dataclass
class EnergyLedger:
    initial_energy: float
    heat_hot: CompensatedCounter = field(default_factory=CompensatedCounter)
    heat_cold: CompensatedCounter = field(default_factory=CompensatedCounter)
    heat_other: CompensatedCounter = field(default_factory=CompensatedCounter)
    work_on: CompensatedCounter = field(default_factory=CompensatedCounter)
    load_output: CompensatedCounter = field(default_factory=CompensatedCounter)
    piston_work_on_gas: CompensatedCounter = field(default_factory=CompensatedCounter)
    support_impulse_x: CompensatedCounter = field(default_factory=CompensatedCounter)
    support_impulse_y: CompensatedCounter = field(default_factory=CompensatedCounter)

    @property
    def heat_into_system(self) -> float:
        return self.heat_hot.value+self.heat_cold.value+self.heat_other.value

    def residual(self, current_energy: float) -> float:
        return current_energy-self.initial_energy-self.heat_into_system-self.work_on.value

    def reservoir_entropy_change(self, t_hot: float, t_cold: float) -> float:
        return -self.heat_hot.value/t_hot-self.heat_cold.value/t_cold
