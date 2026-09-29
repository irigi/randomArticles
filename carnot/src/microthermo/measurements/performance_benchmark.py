"""Reproducible short Carnot solver benchmark; run with python -m."""
from __future__ import annotations

import argparse
import json
import platform
import resource
import time

import numpy as np

from ..api import load_preset
from ..config import RunConfig
from ..core.numeric import (fixed_wall_toi, numba_available, polygon_pair_toi,
                            polygon_wall_gap, polygon_witness_at, reach_mask,
                            wall_lower_bounds)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preset", choices=("carnot_discs", "carnot_triangles"),
                        default="carnot_triangles")
    parser.add_argument("--particles", type=int, default=64)
    parser.add_argument("--duration", type=float, default=.1)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--pair-search", choices=("all", "grid"), default="grid")
    parser.add_argument("--numeric-backend", choices=("auto", "python", "numba"),
                        default="auto")
    parser.add_argument("--wall-search", choices=("bounded", "all"),
                        default="bounded")
    parser.add_argument("--wall-kernel", choices=("auto", "python"),
                        default="auto")
    parser.add_argument("--penetration-kernel", choices=("auto", "python"),
                        default="auto")
    args = parser.parse_args()
    if args.numeric_backend == "numba" and not numba_available():
        parser.error("Numba is not installed; install microthermo[accel]")
    config = RunConfig(preset=args.preset, particles=args.particles, duration=args.duration,
                       seed=args.seed, pair_search=args.pair_search,
                       numeric_backend=args.numeric_backend,
                       wall_search=args.wall_search, wall_kernel=args.wall_kernel,
                       penetration_kernel=args.penetration_kernel)
    actual_backend = ("numba" if args.numeric_backend == "auto" and numba_available()
                      else "python" if args.numeric_backend == "auto"
                      else args.numeric_backend)
    warmup_start = time.perf_counter()
    if actual_backend == "numba":
        reach_mask(np.zeros((2, 2)), np.zeros((2, 2)), np.ones(2), np.zeros(2),
                   np.array([[0, 1]], dtype=np.int64), .01, "numba")
        wall_lower_bounds(np.zeros((2, 2)), np.ones(2),
                          np.zeros((1, 2)), np.array([[1., 0.]]))
        if args.preset == "carnot_triangles":
            triangle = np.array([[0., 1.], [-.8660254, -.5], [.8660254, -.5]])
            zero = np.zeros(2)
            polygon_witness_at(triangle, triangle, zero, np.array([2., 0.]),
                               zero, zero, 0., 0., 0., 0., 0.)
    jit_warmup_seconds = time.perf_counter() - warmup_start
    start = time.perf_counter()
    sim = load_preset(config)
    setup_seconds = time.perf_counter() - start
    if actual_backend == "numba" and args.preset == "carnot_triangles" and args.particles >= 2:
        start = time.perf_counter()
        polygon_pair_toi(sim.world.bodies, 0, 1, min(args.duration, .01), sim.tol, 1024)
        wall = sim.world.walls[0]
        state = sim.world.bodies
        speed = (float(np.linalg.norm(state.vel[0])) +
                 abs(state.omega[0]) * state.radius[0] +
                 wall.speed_bound(0., min(args.duration, .01)))
        fixed_wall_toi(state, 0, wall, wall.point_at(0.), speed,
                       min(args.duration, .01), sim.tol)
        polygon_wall_gap(state, 0, wall.point_at(0.), wall.inward_normal)
        jit_warmup_seconds += time.perf_counter() - start
    start = time.perf_counter()
    sim.advance_to(args.duration)
    solver_seconds = time.perf_counter() - start
    print(json.dumps({
        "preset": args.preset, "particles": args.particles, "seed": args.seed,
        "physical_seconds": args.duration, "setup_seconds": setup_seconds,
        "jit_warmup_seconds": jit_warmup_seconds,
        "solver_seconds": solver_seconds, "events": sim.event_count,
        "events_per_wall_second": sim.event_count / solver_seconds,
        "physical_per_wall_second": args.duration / solver_seconds,
        "peak_process_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "pair_search": args.pair_search,
        "wall_search": args.wall_search,
        "wall_kernel": args.wall_kernel,
        "penetration_kernel": args.penetration_kernel,
        "numeric_backend": actual_backend,
        "particle_radius": sim.world.metadata["particle_radius"],
        "area_fraction_min": sim.world.metadata["particle_area_fraction_min"],
        "max_penetration": sim.max_penetration,
        "ccd_failures": sim.ccd_failures,
        "energy_residual": sim.snapshot().energy_residual,
        "python": platform.python_version(), "numpy": np.__version__,
        "cpu": platform.processor(), "platform": platform.platform(),
    }, indent=2))


if __name__ == "__main__":
    main()
