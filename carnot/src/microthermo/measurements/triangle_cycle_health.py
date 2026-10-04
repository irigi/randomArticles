"""Measure long default-speed triangle runs and optional replay determinism."""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import replace
import gzip
import json
import math
from pathlib import Path
import pickle
import platform
import time

import numpy as np

from ..api import load_preset
from ..config import RunConfig
from ..core.numeric import numba_available
from .performance_benchmark import cpu_model


def save(path: Path, report: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n")
    temporary.replace(path)


def save_branch_checkpoint(path: Path, payload: dict) -> None:
    """Atomically save a trusted local checkpoint after a complete branch."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with gzip.open(temporary, "wb", compresslevel=1) as stream:
        pickle.dump(payload, stream, protocol=pickle.HIGHEST_PROTOCOL)
    temporary.replace(path)


def load_branch_checkpoint(path: Path) -> dict:
    # Pickle is only appropriate for files produced locally by this command.
    with gzip.open(path, "rb") as stream:
        return pickle.load(stream)


def same_result(a, b) -> bool:
    state_a, state_b = a.world.bodies, b.world.bodies
    return (a.events == b.events and a.ledger == b.ledger and
            a.cycle_markers == b.cycle_markers and
            a.failure_diagnostic == b.failure_diagnostic and
            a.ccd_failures == b.ccd_failures and
            all(np.array_equal(getattr(state_a, name), getattr(state_b, name))
                for name in ("pos", "vel", "angle", "omega")))


def compare_same_state(source, compiled, reference, horizon: float) -> dict:
    checkpoint = source.checkpoint()
    compiled.restore(checkpoint)
    reference.restore(checkpoint)
    compiled._pending = None
    reference._pending = None
    left = compiled._earliest(horizon)
    right = reference._earliest(horizon)
    time_delta = (0.0 if left.time == right.time else
                  abs(left.time - right.time)
                  if left.time is not None and right.time is not None and
                  math.isfinite(left.time) and math.isfinite(right.time)
                  else math.inf)
    contact_matches = []
    max_point_delta = 0.0
    max_normal_delta = 0.0
    for a, b in zip(left.contacts, right.contacts):
        max_point_delta = max(max_point_delta, float(np.max(abs(
            np.asarray(a.point) - np.asarray(b.point)))))
        max_normal_delta = max(max_normal_delta, float(np.max(abs(
            np.asarray(a.normal) - np.asarray(b.normal)))))
        contact_matches.append((a.a, a.b, a.boundary, a.feature_a, a.feature_b) ==
                               (b.a, b.b, b.boundary, b.feature_a, b.feature_b))
    matches = (left.status == right.status and left.reason == right.reason and
               time_delta <= 2e-11 and len(left.contacts) == len(right.contacts) and
               all(contact_matches) and max_point_delta <= 2e-10 and
               max_normal_delta <= 2e-10)
    return {"time": compiled.time, "horizon": horizon, "matches": matches,
            "status_numba": left.status, "status_python": right.status,
            "reason_numba": left.reason, "reason_python": right.reason,
            "contacts_numba": len(left.contacts),
            "contacts_python": len(right.contacts),
            "time_delta": time_delta, "max_point_delta": max_point_delta,
            "max_normal_delta": max_normal_delta,
            "contact_features_match": all(contact_matches)}


def measure(config: RunConfig, branches: int, cadence: float | None,
            same_state_queries: bool = False,
            branch_checkpoint: Path | None = None) -> dict:
    duration = branches * math.pi / (2 * config.shaft_speed)
    start = time.perf_counter()
    direct = load_preset(config)
    setup_seconds = time.perf_counter() - start
    branch_rows = []
    query_rows = []
    compiled_probe = load_preset(config) if same_state_queries else None
    reference = (load_preset(replace(config, numeric_backend="python"))
                 if same_state_queries else None)
    solver_seconds = 0.0
    first_branch = 1
    if branch_checkpoint is not None and branch_checkpoint.exists():
        saved = load_branch_checkpoint(branch_checkpoint)
        expected = (1, config.digest(), branches, same_state_queries)
        actual = (saved["version"], saved["config_digest"],
                  saved["branches"], saved["same_state_queries"])
        if actual != expected:
            raise ValueError("branch checkpoint does not match this run")
        first_branch = saved["completed_branch"] + 1
        if not 1 <= first_branch <= branches + 1:
            raise ValueError("invalid completed branch in checkpoint")
        direct.restore(saved["simulation"])
        direct.max_penetration = saved["max_penetration"]
        direct.ccd_refinements = saved["ccd_refinements"]
        direct.ccd_failures = saved["ccd_failures"]
        direct.failure_diagnostic = saved["failure_diagnostic"]
        branch_rows = saved["branch_rows"]
        query_rows = saved["query_rows"]
        setup_seconds = saved["setup_seconds"]
        solver_seconds = saved["solver_seconds"]
        if (len(branch_rows) != first_branch - 1 or
                len(query_rows) != ((first_branch - 1) if same_state_queries else 0) or
                not math.isclose(direct.time,
                                 (first_branch - 1) * math.pi / (2 * config.shaft_speed),
                                 abs_tol=1e-10)):
            raise ValueError("inconsistent branch checkpoint")
        print(f"resumed seed={config.seed} after branch={first_branch - 1}",
              flush=True)
    for branch in range(first_branch, branches + 1):
        if reference is not None:
            start = time.perf_counter()
            direct.advance_to((branch - .5) * duration / branches)
            solver_seconds += time.perf_counter() - start
            query_rows.append(compare_same_state(direct, compiled_probe, reference,
                                                 config.max_horizon))
        start = time.perf_counter()
        direct.advance_to(branch * duration / branches)
        solver_seconds += time.perf_counter() - start
        branch_rows.append({
            "branch": branch,
            "time": direct.time,
            "events": len(direct.events),
            "ccd_failures": direct.ccd_failures,
            "max_penetration": direct.max_penetration,
            "energy_residual": direct.snapshot().energy_residual,
        })
        if branch_checkpoint is not None:
            save_branch_checkpoint(branch_checkpoint, {
                "version": 1, "config_digest": config.digest(),
                "branches": branches, "same_state_queries": same_state_queries,
                "completed_branch": branch, "simulation": direct.checkpoint(),
                "max_penetration": direct.max_penetration,
                "ccd_refinements": direct.ccd_refinements,
                "ccd_failures": direct.ccd_failures,
                "failure_diagnostic": direct.failure_diagnostic,
                "branch_rows": branch_rows, "query_rows": query_rows,
                "setup_seconds": setup_seconds,
                "solver_seconds": solver_seconds,
            })
            print(f"seed={config.seed} branch={branch}/{branches} "
                  f"events={len(direct.events)}", flush=True)
    row = {
        "seed": config.seed,
        "config_digest": config.digest(),
        "setup_seconds": setup_seconds,
        "solver_seconds": solver_seconds,
        "physical_per_wall_second": duration / solver_seconds,
        "duration": duration,
        "events": len(direct.events),
        "event_kinds": dict(Counter(event.kind for event in direct.events)),
        "completed_cycles": len(direct.cycle_markers) - 1,
        "ccd_refinements": direct.ccd_refinements,
        "ccd_failures": direct.ccd_failures,
        "max_penetration": direct.max_penetration,
        "energy_residual": direct.snapshot().energy_residual,
        "ledger": {name: getattr(direct.ledger, name).value for name in (
            "heat_hot", "heat_cold", "heat_other", "work_on",
            "load_output", "piston_work_on_gas")},
        "particle_area_fraction_min": direct.world.metadata[
            "particle_area_fraction_min"],
        "radius": float(direct.world.bodies.radius[0]),
        "branches": branch_rows,
    }
    if reference is not None:
        row["same_state_queries"] = query_rows
    if cadence is not None:
        sampled = load_preset(config)
        checkpoint = None
        start = time.perf_counter()
        for target in np.arange(cadence, duration, cadence):
            sampled.advance_to(float(target))
            if checkpoint is None and target > duration / 2:
                checkpoint = sampled.checkpoint()
        sampled.advance_to(duration)
        row["sampling"] = {
            "cadence": cadence,
            "sampled_seconds": time.perf_counter() - start,
            "matches_direct": same_result(sampled, direct),
        }
        if checkpoint is not None:
            sampled.restore(checkpoint)
            start = time.perf_counter()
            sampled.advance_to(duration)
            row["sampling"]["replay_seconds"] = time.perf_counter() - start
            row["sampling"]["replay_matches_direct"] = same_result(
                sampled, direct)
    return row


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--particles", type=int, default=96)
    parser.add_argument("--seeds", type=int, nargs="+", default=[123, 124, 125])
    parser.add_argument("--cycles", type=int, default=2)
    parser.add_argument("--branches", type=int,
                        help="number of Carnot branches; overrides --cycles")
    parser.add_argument("--shaft-speed", type=float, default=.15)
    parser.add_argument("--verify-seed", type=int, default=123,
                        help="seed to rerun with sampling and checkpoint replay")
    parser.add_argument("--no-verify", action="store_true",
                        help="skip the sampled run and checkpoint replay")
    parser.add_argument("--same-state-queries", action="store_true",
                        help="compare Python and Numba next collisions at branch midpoints")
    parser.add_argument("--sample-cadence", type=float, default=.83)
    parser.add_argument("--branch-checkpoint", type=Path,
                        help="save and resume trusted local branch-end state; one seed only")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if (args.particles < 1 or args.cycles < 1 or
            (args.branches is not None and args.branches < 1) or
            args.shaft_speed <= 0 or
            args.sample_cadence <= 0 or not args.seeds or
            len(set(args.seeds)) != len(args.seeds) or
            (not args.no_verify and args.verify_seed not in args.seeds) or
            (args.branch_checkpoint is not None and len(args.seeds) != 1)):
        parser.error("require positive settings, unique seeds, and verify-seed in seeds")
    if not numba_available():
        parser.error("Numba is required for this long-run measurement")
    branches = args.branches if args.branches is not None else 4 * args.cycles
    report = {
        "schema_version": 1,
        "complete": False,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "cpu": cpu_model(),
        "settings": {"particles": args.particles, "seeds": args.seeds,
                     "cycles": args.cycles if args.branches is None else None,
                     "branches": branches, "shaft_speed": args.shaft_speed,
                     "verify_seed": None if args.no_verify else args.verify_seed,
                     "sample_cadence": args.sample_cadence,
                     "same_state_queries": args.same_state_queries,
                     "numeric_backend": "numba", "pair_search": "sweep"},
        "runs": [],
    }
    save(args.output, report)
    for seed in args.seeds:
        config = RunConfig(preset="carnot_triangles", particles=args.particles,
                           seed=seed, shaft_speed=args.shaft_speed,
                           pair_search="sweep", numeric_backend="numba")
        try:
            row = measure(config, branches,
                          args.sample_cadence if (not args.no_verify and
                                                  seed == args.verify_seed)
                          else None, args.same_state_queries,
                          args.branch_checkpoint)
        except Exception as exc:
            report["runs"].append({"seed": seed, "error": str(exc),
                                   "error_type": type(exc).__name__})
            save(args.output, report)
            raise
        report["runs"].append(row)
        save(args.output, report)
        print(f"seed={seed} events={row['events']} "
              f"failures={row['ccd_failures']} "
              f"penetration={row['max_penetration']:.3g} "
              f"residual={row['energy_residual']:.3g}", flush=True)
    report["complete"] = True
    save(args.output, report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
