"""Throughput of the edmd kernel against the reference engine; run with python -m.

Both engines run the same world: smooth discs at a fixed area fraction in a
box whose left wall is thermal, like `gas_box`, scaled with N so the density
stays constant. Each engine runs for a wall-clock budget after a warm-up
that excludes Numba compilation.
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import time
from pathlib import Path

import numpy as np

from ..core.boundaries import BoundaryKind
from ..core.state import BodyArrays
from ..experiments.common import box_walls, particle_specs
from ..runner.edmd_simulation import EdmdSimulation
from ..runner.simulation import Simulation, World
from .performance_benchmark import cpu_model


def scaled_world(particles: int, area_fraction: float, radius: float, seed: int) -> World:
    area = particles*math.pi*radius**2/area_fraction
    height = math.sqrt(area/2)
    width = 2*height
    rng = np.random.default_rng(seed)
    specs = particle_specs(particles, width, height, radius, 1., rng)
    return World(BodyArrays.from_specs(specs),
                 box_walls(width, height, BoundaryKind.HOT, 1.),
                 metadata={"name": "scaled_gas_box", "width": width, "height": height})


def measure(make, budget: float, step: float) -> dict:
    sim = make()
    sim.advance_to(step)                       # warm-up, excluded
    t0, e0 = sim.time, sim.event_count
    start = time.perf_counter()
    while time.perf_counter()-start < budget:
        sim.advance_to(sim.time+step)
        if hasattr(sim, "drain_events"):
            sim.drain_events()                 # as precalculation does
    wall = time.perf_counter()-start
    events = sim.event_count-e0
    physical = sim.time-t0
    residual = sim.snapshot().energy_residual
    return {"wall_seconds": wall, "events": events, "physical_seconds": physical,
            "events_per_second": events/wall, "physical_per_wall_second": physical/wall,
            "energy_residual": residual}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--particles", type=int, nargs="+", default=[200, 1000, 5000])
    parser.add_argument("--area-fraction", type=float, default=.05)
    parser.add_argument("--radius", type=float, default=.02)
    parser.add_argument("--budget", type=float, default=20.,
                        help="wall-clock seconds per engine and size")
    parser.add_argument("--reference-max", type=int, default=1000,
                        help="skip the reference engine above this many particles")
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--output")
    args = parser.parse_args()
    rows = []
    for n in args.particles:
        world = lambda: scaled_world(n, args.area_fraction, args.radius, args.seed)
        # (factory, sample step): .05 is the usual sampling cadence; the
        # long step shows the kernel itself, without per-sample overhead.
        runs = {
            "edmd": (lambda: EdmdSimulation(world(), args.seed), .05),
            "edmd_no_event_log": (lambda: EdmdSimulation(world(), args.seed,
                                                         record_events=False), .05),
            "edmd_kernel_only": (lambda: EdmdSimulation(world(), args.seed,
                                                        record_events=False), 2.),
        }
        if n <= args.reference_max:
            runs["reference"] = (lambda: Simulation(world(), args.seed, .05,
                                                    pair_search="grid"), .05)
        for engine, (make, step) in runs.items():
            row = {"particles": n, "engine": engine, "sample_step": step,
                   **measure(make, args.budget, step)}
            rows.append(row)
            print(json.dumps(row), flush=True)
    for n in args.particles:
        by = {r["engine"]: r for r in rows if r["particles"] == n}
        if "reference" in by:
            for engine in ("edmd", "edmd_no_event_log", "edmd_kernel_only"):
                by[engine]["speedup_vs_reference"] = (by[engine]["events_per_second"] /
                                                      by["reference"]["events_per_second"])
    report = {"date": time.strftime("%Y-%m-%d"), "platform": platform.platform(),
              "cpu": cpu_model(), "python": platform.python_version(),
              "numpy": np.__version__, "area_fraction": args.area_fraction,
              "radius": args.radius, "budget_seconds": args.budget, "results": rows}
    text = json.dumps(report, indent=2)+"\n"
    if args.output:
        Path(args.output).write_text(text)
    print(text)


if __name__ == "__main__":
    main()
