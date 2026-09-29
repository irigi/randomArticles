"""Reproducible controlled-shaft speed/seed measurements.

This is a descriptive data collection tool. It makes no convergence claim.
"""

from __future__ import annotations

from dataclasses import asdict
import math
import statistics

from ..api import load_preset
from ..config import RunConfig
from ..runner.simulation import NumericalFailure
from .cycles import assess_efficiency, summarize_cycle
from .ideal import ideal_carnot_reference
from .stationarity import cycle_drift, gas_energy_drift


BRANCH_NAMES=("hot","adiabatic_expansion","cold","adiabatic_compression")


def _branch_record(name, cycle_index, start, end, start_gas, end_gas, events):
    thermal=[event for event in events if event.kind in ("hot","cold")]
    heats=[event.heat_into_system for event in thermal]
    mode_in=[event.metadata["incoming_normal_mode_energy"] for event in thermal]
    mode_out=[event.metadata["outgoing_normal_mode_energy"] for event in thermal]
    return {"cycle":cycle_index+1,"branch":name,
            "start_time":start.time,"end_time":end.time,
            "temperature_trans_start":start.translational_temperature,
            "temperature_trans_end":end.translational_temperature,
            "temperature_rot_start":start.rotational_temperature,
            "temperature_rot_end":end.rotational_temperature,
            "gas_energy_start":start_gas,
            "gas_energy_end":end_gas,
            "thermal_contacts":len(thermal),
            "heat_into_gas":math.fsum(heats),
            "heat_added":math.fsum(value for value in heats if value>0),
            "heat_removed":math.fsum(value for value in heats if value<0),
            "normal_mode_energy_in_total":math.fsum(mode_in),
            "normal_mode_energy_out_total":math.fsum(mode_out),
            "normal_mode_heat_reconciliation":(
                math.fsum(mode_out)-math.fsum(mode_in)-math.fsum(heats))}


def run_speed_study(*, preset: str, speeds: tuple[float, ...],
                    seeds: tuple[int, ...], cycles: int, particles: int,
                    max_horizon: float = .05, transient_cycles: int = 2,
                    efficiency_min_cycles: int = 8) -> dict:
    """Measure complete controlled cycles at each speed and independent seed."""
    if preset not in ("carnot_discs", "carnot_triangles"):
        raise ValueError("speed study requires a Carnot preset")
    if not speeds or not seeds or len(set(speeds)) != len(speeds) or len(set(seeds)) != len(seeds):
        raise ValueError("require nonempty, unique speeds and seeds")
    if cycles <= 0 or particles <= 0 or max_horizon <= 0:
        raise ValueError("cycles, particles, and max_horizon must be positive")
    if any(not math.isfinite(speed) or speed <= 0 for speed in speeds):
        raise ValueError("shaft speeds must be finite and positive")
    if transient_cycles < 0 or efficiency_min_cycles < 4:
        raise ValueError("invalid efficiency cycle gate")

    groups=[]
    for speed in speeds:
        runs=[]
        ideal_efficiency=None
        ideal_work=None
        for seed in seeds:
            duration=cycles*2*math.pi/speed
            config=RunConfig(preset=preset,seed=seed,particles=particles,
                             duration=duration,cycles=cycles,max_horizon=max_horizon,
                             shaft_mode="controlled",shaft_speed=speed,
                             transient_cycles=transient_cycles,
                             efficiency_min_cycles=efficiency_min_cycles)
            simulation=load_preset(config)
            reference=ideal_carnot_reference(simulation)
            if reference is not None:
                ideal_efficiency=reference.efficiency
                ideal_work=reference.work_by_gas
            try:
                cam=simulation.world.mechanism.cam
                branch_records=[]
                start_snapshot=simulation.snapshot()
                start_gas=simulation.world.bodies.kinetic_energy()
                event_cursor=0
                for cycle_index in range(cycles):
                    for branch_index,name in enumerate(BRANCH_NAMES):
                        target=(cycle_index*2*math.pi+
                                float(cam.boundaries[branch_index+1]))/speed
                        end_snapshot=simulation.advance_to(target)
                        end_gas=simulation.world.bodies.kinetic_energy()
                        branch_records.append(_branch_record(
                            name,cycle_index,start_snapshot,end_snapshot,
                            start_gas,end_gas,
                            simulation.events[event_cursor:]))
                        start_snapshot=end_snapshot
                        start_gas=end_gas
                        event_cursor=len(simulation.events)
            except NumericalFailure as exc:
                runs.append({"seed":seed,"config_hash":config.digest(),
                             "status":"numerical_failure","error":str(exc),
                             "completed_cycles":max(0,len(simulation.cycle_markers)-1)})
                continue
            summaries=[summarize_cycle(start,end) for start,end in zip(
                simulation.cycle_markers[:-1],simulation.cycle_markers[1:])]
            report=assess_efficiency(summaries,transient_cycles,
                                     efficiency_min_cycles)
            eligible=summaries[transient_cycles:]
            eligible_branches=[record for record in branch_records
                               if record["cycle"]>transient_cycles]
            hot_branches=[record for record in eligible_branches
                          if record["branch"]=="hot"]
            hot_heat=math.fsum(x.heat_hot for x in eligible)
            hot_branch_heat=math.fsum(x["heat_into_gas"] for x in hot_branches)
            hot_contacts=sum(x["thermal_contacts"] for x in hot_branches)
            hot_mode_in=math.fsum(x["normal_mode_energy_in_total"]
                                 for x in hot_branches)
            hot_mode_out=math.fsum(x["normal_mode_energy_out_total"]
                                  for x in hot_branches)
            runs.append({
                "seed":seed,"config_hash":config.digest(),
                "status":report.status,"completed_cycles":len(summaries),
                "event_count":simulation.event_count,
                "max_penetration":simulation.max_penetration,
                "ccd_failures":simulation.ccd_failures,
                "first_law_residual":simulation.snapshot().energy_residual,
                "post_transient_hot_heat":hot_heat,
                "hot_branch_heat_reconciliation":hot_branch_heat-hot_heat,
                "hot_branch_contacts":hot_contacts,
                "mean_incoming_hot_normal_mode_energy":(
                    hot_mode_in/hot_contacts if hot_contacts else None),
                "mean_outgoing_hot_normal_mode_energy":(
                    hot_mode_out/hot_contacts if hot_contacts else None),
                "hot_normal_mode_heat_reconciliation":(
                    hot_mode_out-hot_mode_in-hot_heat),
                "hot_branch_heat_added":math.fsum(x["heat_added"]
                                                   for x in hot_branches),
                "hot_branch_heat_removed":math.fsum(x["heat_removed"]
                                                     for x in hot_branches),
                "mean_pre_hot_temperature":(statistics.mean(
                    x["temperature_trans_start"] for x in hot_branches)
                    if hot_branches else None),
                "mean_post_hot_temperature":(statistics.mean(
                    x["temperature_trans_end"] for x in hot_branches)
                    if hot_branches else None),
                "post_transient_external_output":math.fsum(
                    x.external_output for x in eligible),
                "post_transient_storage_change":math.fsum(
                    x.delta_gas_energy+x.delta_apparatus_energy for x in eligible),
                "efficiency_report":asdict(report),
                "gas_energy_drift":gas_energy_drift(
                    [marker.gas_energy for marker in
                     simulation.cycle_markers[1:]],transient_cycles),
                "hot_heat_drift":cycle_drift(
                    [cycle.heat_hot for cycle in summaries],transient_cycles),
                "branch_diagnostics":branch_records,
                "cycles": [asdict(x) for x in summaries],
            })
        values=[run["efficiency_report"]["estimate"]["value"] for run in runs
                if run["status"] == "ready"]
        groups.append({
            "shaft_speed":speed,"nominal_cycle_duration":2*math.pi/speed,
            "ideal_efficiency":ideal_efficiency,
            "ideal_work_by_gas_per_cycle":ideal_work,
            "ready_seed_count":len(values),
            "mean_measured_efficiency":(statistics.mean(values)
                                        if len(values)>=2 else None),
            "seed_standard_error":(statistics.stdev(values)/math.sqrt(len(values))
                                   if len(values)>=2 else None),
            "runs":runs,
        })
    return {"preset":preset,"shaft_mode":"controlled", "particles":particles,
            "requested_cycles":cycles,"seeds":list(seeds),
            "transient_cycles":transient_cycles,
            "efficiency_min_cycles":efficiency_min_cycles,
            "uncertainty_note":"Seed standard error is descriptive; no convergence test is applied.",
            "speeds":groups}
