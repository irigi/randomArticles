"""Reproducible short Carnot solver benchmark; run with python -m."""
from __future__ import annotations

import argparse
import json
import platform
import resource
import time
from pathlib import Path

import numpy as np

from ..api import load_preset
from ..config import RunConfig
from ..core.broadphase import swept_pairs_sweep
from ..core.numeric import (controlled_cam_disc_toi,
                            controlled_cam_triangle_toi, disc_pairs_toi,
                            disc_penetration,
                            fixed_disc_walls_toi,
                            fixed_wall_toi, numba_available,
                            polygon_pair_toi,
                            polygon_pairs_toi,
                            polygon_wall_gap, polygon_witness_at, reach_mask,
                            swept_circle_mask,
                            triangle_pair_gaps, triangle_penetration,
                            wall_lower_bounds, wall_reach_mask)
from ..core.state import BodyArrays, BodySpec, Shape


def cpu_model() -> str:
    info = Path("/proc/cpuinfo")
    if info.exists():
        for line in info.read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    return platform.processor() or platform.machine()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preset", choices=("carnot_discs", "carnot_triangles"),
                        default="carnot_triangles")
    parser.add_argument("--particles", type=int, default=64)
    parser.add_argument("--duration", type=float, default=.1)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--shaft-speed", type=float, default=.15)
    parser.add_argument("--pair-search", choices=("all", "grid", "sweep"), default="grid")
    parser.add_argument("--numeric-backend", choices=("auto", "python", "numba"),
                        default="auto")
    parser.add_argument("--wall-search", choices=("bounded", "all"),
                        default="bounded")
    parser.add_argument("--wall-kernel", choices=("auto", "python"),
                        default="auto")
    parser.add_argument("--penetration-kernel", choices=("auto", "python"),
                        default="auto")
    parser.add_argument("--pair-kernel", choices=("auto", "scalar"), default="auto")
    parser.add_argument("--cam-kernel", choices=("auto", "python"), default="auto")
    args = parser.parse_args()
    if args.numeric_backend == "numba" and not numba_available():
        parser.error("Numba is not installed; install microthermo[accel]")
    config = RunConfig(preset=args.preset, particles=args.particles, duration=args.duration,
                       seed=args.seed, shaft_speed=args.shaft_speed,
                       pair_search=args.pair_search,
                       numeric_backend=args.numeric_backend,
                       wall_search=args.wall_search, wall_kernel=args.wall_kernel,
                       penetration_kernel=args.penetration_kernel,
                       pair_kernel=args.pair_kernel, cam_kernel=args.cam_kernel)
    actual_backend = ("numba" if args.numeric_backend == "auto" and numba_available()
                      else "python" if args.numeric_backend == "auto"
                      else args.numeric_backend)
    warmup_start = time.perf_counter()
    if actual_backend == "numba":
        if (args.preset == "carnot_discs" and
                args.penetration_kernel == "auto" and
                args.wall_search == "bounded"):
            disc_penetration(np.zeros((1, 2)), np.ones(1),
                             np.empty((0, 2), dtype=np.int64),
                             np.empty((0, 2)), np.empty((0, 2)))
        if args.pair_search == "sweep":
            swept_pairs_sweep(np.zeros((2, 2)), np.zeros((2, 2)),
                              np.ones(2), .01, "numba")
        reach_mask(np.zeros((2, 2)), np.zeros((2, 2)), np.ones(2), np.zeros(2),
                   np.array([[0, 1]], dtype=np.int64), .01, "numba")
        if args.pair_kernel == "auto":
            swept_circle_mask(np.zeros((2, 2)), np.zeros((2, 2)), np.ones(2),
                              np.array([[0, 1]], dtype=np.int64), .01, 1e-10,
                              "numba")
        wall_lower_bounds(np.zeros((2, 2)), np.ones(2),
                          np.zeros((1, 2)), np.array([[1., 0.]]))
        wall_reach_mask(np.zeros((2, 2)), np.zeros((2, 2)), np.ones(2),
                        np.zeros((1, 2)), np.array([[1., 0.]]), np.ones(1),
                        np.zeros(1, dtype=np.bool_), .01, 1e-10, "numba")
        if args.preset == "carnot_triangles":
            triangle = np.array([[0., 1.], [-.8660254, -.5], [.8660254, -.5]])
            zero = np.zeros(2)
            polygon_witness_at(triangle, triangle, zero, np.array([2., 0.]),
                               zero, zero, 0., 0., 0., 0., 0.)
            if args.penetration_kernel == "auto":
                state = BodyArrays.from_specs([
                    BodySpec((0., 0.), (0., 0.), shape=Shape.TRIANGLE),
                    BodySpec((1., 0.), (0., 0.), shape=Shape.TRIANGLE)])
                triangle_pair_gaps(state, np.array([[0, 1]], dtype=np.int64))
                triangle_penetration(
                    state, np.array([[0, 1]], dtype=np.int64),
                    np.array([[0., -2.]]), np.array([[0., 1.]]))
    jit_warmup_seconds = time.perf_counter() - warmup_start
    start = time.perf_counter()
    sim = load_preset(config)
    setup_seconds = time.perf_counter() - start
    if actual_backend == "numba" and args.preset == "carnot_triangles" and args.particles >= 2:
        start = time.perf_counter()
        polygon_pair_toi(sim.world.bodies, 0, 1, min(args.duration, .01), sim.tol, 1024)
        if args.pair_kernel == "auto":
            polygon_pairs_toi(sim.world.bodies, np.array([[0, 1]], dtype=np.int64),
                              min(args.duration, .01), sim.tol, 1024)
        wall = sim.world.walls[0]
        state = sim.world.bodies
        speed = (float(np.linalg.norm(state.vel[0])) +
                 abs(state.omega[0]) * state.radius[0] +
                 wall.speed_bound(0., min(args.duration, .01)))
        fixed_wall_toi(state, 0, wall, wall.point_at(0.), speed,
                       min(args.duration, .01), sim.tol)
        polygon_wall_gap(state, 0, wall.point_at(0.), wall.inward_normal)
        jit_warmup_seconds += time.perf_counter() - start
    if (actual_backend == "numba" and args.preset == "carnot_discs" and
            args.pair_kernel == "auto" and args.particles >= 2):
        start = time.perf_counter()
        disc_pairs_toi(sim.world.bodies, np.array([[0, 1]], dtype=np.int64),
                       min(args.duration, .01), sim.tol)
        jit_warmup_seconds += time.perf_counter() - start
    if (actual_backend == "numba" and args.preset == "carnot_discs" and
            args.cam_kernel == "auto"):
        start = time.perf_counter()
        wall = sim.world.walls[3]
        controlled_cam_disc_toi(sim.world.bodies,
                                np.zeros(args.particles, dtype=np.bool_),
                                sim.world.mechanism,
                                wall.speed_bound(0., min(args.duration, .01)),
                                min(args.duration, .01), sim.tol)
        jit_warmup_seconds += time.perf_counter() - start
    if (actual_backend == "numba" and args.preset == "carnot_triangles" and
            args.cam_kernel == "auto"):
        start = time.perf_counter()
        wall = sim.world.walls[3]
        controlled_cam_triangle_toi(
            sim.world.bodies, np.zeros(args.particles, dtype=np.bool_),
            sim.world.mechanism,
            wall.speed_bound(0., min(args.duration, .01)),
            min(args.duration, .01), sim.tol)
        jit_warmup_seconds += time.perf_counter() - start
    if (actual_backend == "numba" and args.preset == "carnot_discs" and
            args.wall_kernel == "auto" and args.wall_search == "bounded"):
        start = time.perf_counter()
        walls = sim.world.walls[:3]
        fixed_disc_walls_toi(
            sim.world.bodies, np.zeros((args.particles, len(walls)), dtype=np.bool_),
            np.asarray([wall.point_at(0.) for wall in walls]),
            np.asarray([wall.inward_normal for wall in walls]),
            np.asarray([wall.velocity for wall in walls]),
            np.asarray([wall.speed_bound(0., .01) for wall in walls]),
            .01, sim.tol)
        jit_warmup_seconds += time.perf_counter() - start
    if actual_backend == "numba":
        # Warm the complete event search for both shapes. Individual kernel
        # calls above cannot cover every scheduler specialization.
        start = time.perf_counter()
        sim._earliest(sim.max_horizon)
        jit_warmup_seconds += time.perf_counter() - start
    start = time.perf_counter()
    sim.advance_to(args.duration)
    solver_seconds = time.perf_counter() - start
    print(json.dumps({
        "preset": args.preset, "particles": args.particles, "seed": args.seed,
        "shaft_speed": args.shaft_speed,
        "physical_seconds": args.duration, "setup_seconds": setup_seconds,
        "jit_warmup_seconds": jit_warmup_seconds,
        "startup_seconds": setup_seconds + jit_warmup_seconds,
        "solver_seconds": solver_seconds, "events": sim.event_count,
        "events_per_wall_second": sim.event_count / solver_seconds,
        "physical_per_wall_second": args.duration / solver_seconds,
        "wall_seconds_per_physical_second": solver_seconds / args.duration,
        "peak_process_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "pair_search": args.pair_search,
        "wall_search": args.wall_search,
        "wall_kernel": args.wall_kernel,
        "penetration_kernel": args.penetration_kernel,
        "pair_kernel": args.pair_kernel,
        "cam_kernel": args.cam_kernel,
        "numeric_backend": actual_backend,
        "particle_radius": sim.world.metadata["particle_radius"],
        "area_fraction_min": sim.world.metadata["particle_area_fraction_min"],
        "max_penetration": sim.max_penetration,
        "ccd_failures": sim.ccd_failures,
        "energy_residual": sim.snapshot().energy_residual,
        "python": platform.python_version(), "numpy": np.__version__,
        "cpu": cpu_model(), "platform": platform.platform(),
    }, indent=2))


if __name__ == "__main__":
    main()
