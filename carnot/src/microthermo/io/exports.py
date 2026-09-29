from __future__ import annotations

from dataclasses import asdict
import csv
import json
from pathlib import Path
from typing import Any
import numpy as np
from ..measurements.cycles import assess_efficiency, summarize_cycle
from ..measurements.pressure_area import compare_pressure_area


def _json(value: Any):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if hasattr(value, "value"): return value.value
    raise TypeError(type(value).__name__)


def _finite(value: float | None):
    return value if value is None or np.isfinite(value) else None


def export_run(output: str | Path, config, simulation) -> Path:
    path = Path(output)
    path.mkdir(parents=True, exist_ok=True)
    resolved = asdict(config)
    resolved["config_hash"] = config.digest()
    (path/"config.json").write_text(json.dumps(resolved, indent=2, default=_json)+"\n")
    with (path/"events.jsonl").open("w") as f:
        for event in simulation.events:
            f.write(json.dumps(event.as_dict(), default=_json)+"\n")
    with (path/"samples.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time", "energy", "temperature_trans", "temperature_rot",
                         "branch", "events", "energy_residual", "occupancy"])
        for s in simulation.samples:
            writer.writerow([s.time, s.energy, s.translational_temperature,
                             s.rotational_temperature, s.branch or "", s.event_count,
                             s.energy_residual, ";".join(map(str, s.occupancy))])
    final = simulation.snapshot()
    cycles = []
    cycle_summaries=[]
    for start,end in zip(simulation.cycle_markers[:-1],simulation.cycle_markers[1:]):
        cycle=summarize_cycle(start,end)
        cycle_summaries.append(cycle)
        record = asdict(cycle)
        comparison = compare_pressure_area(simulation,start,end)
        record["pressure_area"] = asdict(comparison) if comparison else None
        cycles.append(record)
    (path/"cycles.json").write_text(json.dumps(cycles,indent=2,default=_json)+"\n")
    efficiency_report=assess_efficiency(
        cycle_summaries,config.transient_cycles,config.efficiency_min_cycles)
    summary = {
        "preset": config.preset, "seed": config.seed, "physical_time": simulation.time,
        "events": simulation.event_count, "energy": final.energy,
        "energy_residual": final.energy_residual,
        "temperature_trans": _finite(final.translational_temperature),
        "temperature_rot": _finite(final.rotational_temperature),
        "heat_hot": simulation.ledger.heat_hot.value,
        "heat_cold": simulation.ledger.heat_cold.value,
        "work_on": simulation.ledger.work_on.value,
        "load_output": simulation.ledger.load_output.value,
        "completed_cycles": len(cycles),
        "last_cycle": cycles[-1] if cycles else None,
        "efficiency": (efficiency_report.estimate.value
                       if efficiency_report.estimate else None),
        "efficiency_report": asdict(efficiency_report),
        "occupancy": final.occupancy, "ccd_refinements": simulation.ccd_refinements,
        "ccd_failures": simulation.ccd_failures,
        "rng": simulation.rng.bit_generator.__class__.__name__,
        "metadata": {k:v for k,v in simulation.world.metadata.items() if k != "clock"},
    }
    (path/"summary.json").write_text(json.dumps(summary, indent=2, default=_json)+"\n")
    return path
