"""Exact-boundary Carnot cycle accounting."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class CycleMarker:
    time: float
    phase: float
    total_energy: float
    gas_energy: float
    apparatus_energy: float
    heat_hot: float
    heat_cold: float
    heat_other: float
    work_on: float
    piston_work_on_gas: float
    motor_work: float
    load_output: float
    event_count: int


@dataclass(frozen=True)
class CycleSummary:
    start_time: float
    end_time: float
    start_phase: float
    end_phase: float
    heat_hot: float
    heat_cold: float
    heat_other: float
    piston_work_on_gas: float
    motor_work: float
    load_output: float
    delta_gas_energy: float
    delta_apparatus_energy: float
    total_first_law_residual: float
    gas_first_law_residual: float
    event_count: int

    @property
    def external_output(self) -> float:
        return self.load_output-self.motor_work


def summarize_cycle(start: CycleMarker, end: CycleMarker) -> CycleSummary:
    hot=end.heat_hot-start.heat_hot
    cold=end.heat_cold-start.heat_cold
    other=end.heat_other-start.heat_other
    gas_work=end.piston_work_on_gas-start.piston_work_on_gas
    delta_gas=end.gas_energy-start.gas_energy
    delta_apparatus=end.apparatus_energy-start.apparatus_energy
    return CycleSummary(
        start.time,end.time,start.phase,end.phase,hot,cold,other,gas_work,
        end.motor_work-start.motor_work,end.load_output-start.load_output,
        delta_gas,delta_apparatus,
        end.total_energy-start.total_energy-(hot+cold+other)-
        (end.work_on-start.work_on),
        delta_gas-(hot+cold+other)-gas_work,
        end.event_count-start.event_count)


@dataclass(frozen=True)
class EfficiencyEstimate:
    """Ratio of summed net external work to summed hot heat."""

    value: float
    block_standard_error: float
    sample_count: int
    block_count: int
    direction: int
    hot_heat: float
    external_output: float
    storage_change: float


@dataclass(frozen=True)
class EfficiencyReport:
    status: str
    completed_cycles: int
    matched_cycles: int
    required_cycles: int
    transient_cycles: int
    estimate: EfficiencyEstimate | None


def assess_efficiency(cycles: list[CycleSummary] | tuple[CycleSummary, ...],
                      transient_cycles: int = 2,
                      min_cycles: int = 8) -> EfficiencyReport:
    """Gate a descriptive efficiency on consecutive full cycles and hot input.

    The uncertainty is a four-block delete-one-block jackknife standard error
    on contiguous cycle blocks. It describes between-block variability, not a
    confidence interval or proof of quasistatic behavior.
    """
    if transient_cycles < 0 or min_cycles < 4:
        raise ValueError("require nonnegative transients and at least four cycles")
    matched: list[CycleSummary] = []
    direction = 0
    for cycle in cycles:
        span = cycle.end_phase-cycle.start_phase
        sign = (1 if math.isclose(span,2*math.pi,abs_tol=1e-8) else
                -1 if math.isclose(span,-2*math.pi,abs_tol=1e-8) else 0)
        if not sign:
            matched=[]
            direction=0
            continue
        if (sign != direction or not matched or
            not math.isclose(cycle.start_phase,matched[-1].end_phase,abs_tol=1e-8) or
            not math.isclose(cycle.start_time,matched[-1].end_time,abs_tol=1e-8)):
            matched=[]
        matched.append(cycle)
        direction=sign
    eligible=matched[transient_cycles:]
    base=dict(completed_cycles=len(cycles),matched_cycles=len(eligible),
              required_cycles=min_cycles,transient_cycles=transient_cycles)
    if len(eligible) < min_cycles:
        return EfficiencyReport("awaiting_matched_cycles",**base,estimate=None)
    hot=[cycle.heat_hot for cycle in eligible]
    work=[cycle.external_output for cycle in eligible]
    total_hot=math.fsum(hot)
    total_work=math.fsum(work)
    heat_floor=max(1e-12,1e-10*math.fsum(abs(q) for q in hot))
    if total_hot <= heat_floor:
        return EfficiencyReport("insufficient_hot_heat",**base,estimate=None)
    blocks=4
    leave_out=[]
    for block in range(blocks):
        first=(block*len(eligible))//blocks
        end=((block+1)*len(eligible))//blocks
        block_hot=math.fsum(hot[first:end])
        if block_hot <= heat_floor:
            return EfficiencyReport("insufficient_hot_heat",**base,estimate=None)
        remaining_hot=math.fsum(hot[:first]+hot[end:])
        if remaining_hot <= heat_floor:
            return EfficiencyReport("insufficient_hot_heat",**base,estimate=None)
        remaining_work=math.fsum(work[:first]+work[end:])
        leave_out.append(remaining_work/remaining_hot)
    mean=math.fsum(leave_out)/blocks
    standard_error=math.sqrt((blocks-1)/blocks*
        math.fsum((value-mean)**2 for value in leave_out))
    storage=math.fsum(cycle.delta_gas_energy+cycle.delta_apparatus_energy
                      for cycle in eligible)
    estimate=EfficiencyEstimate(total_work/total_hot,standard_error,
        len(eligible),blocks,direction,total_hot,total_work,storage)
    return EfficiencyReport("ready",**base,estimate=estimate)
