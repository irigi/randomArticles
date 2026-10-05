"""Print the per-cycle Carnot ledger for one triangle configuration.

Usage (from the repository root):
    .venv/bin/python examples/carnot_cycle_pilot.py NAME PARTICLES CYCLES '{"shaft_speed": 0.0375,
        "hot_jacket": true, "cold_jacket": true, "initial_temperature": 1.5}'

The JSON holds RunConfig overrides; shaft_speed defaults to 0.15.
W_by_gas is work delivered by the gas (positive for an engine).
"""
import json
import math
import sys
import time

from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.measurements.cycles import summarize_cycle

name, particles, cycles = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
extra = json.loads(sys.argv[4]) if len(sys.argv) > 4 else {}
speed = extra.pop("shaft_speed", .15)
period = 2*math.pi/speed
config = RunConfig(preset="carnot_triangles", particles=particles, seed=123,
                   shaft_speed=speed, duration=cycles*period, max_horizon=.05,
                   **extra)
simulation = load_preset(config)
start = time.perf_counter()
for k in range(1, cycles+1):
    simulation.advance_to(k*period)
    simulation.drain_events()
    markers = simulation.cycle_markers
    cycle = summarize_cycle(markers[-2], markers[-1])
    work = -cycle.piston_work_on_gas
    print(json.dumps({
        "name": name, "cycle": k, "QH": round(cycle.heat_hot, 2),
        "QC": round(cycle.heat_cold, 2), "W_by_gas": round(work, 2),
        "eta": round(work/cycle.heat_hot, 3) if cycle.heat_hot > 0 else None,
        "T_end": round(simulation.snapshot().translational_temperature, 3),
        "wall_s": round(time.perf_counter()-start, 1)}), flush=True)
