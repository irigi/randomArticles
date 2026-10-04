"""Run the Carnot throughput matrix in isolated processes and save JSON."""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from pathlib import Path


STANDARD_DURATIONS = {16: .3, 32: .2, 64: .1, 128: .04, 200: .02, 500: .005}
PRESETS = ("carnot_discs", "carnot_triangles")


def parse_csv(value: str, *, positive: bool = False) -> list[int]:
    try:
        values = [int(part.strip()) for part in value.split(",")]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from exc
    if not values or any(value < (1 if positive else 0) for value in values):
        raise argparse.ArgumentTypeError("expected positive integers")
    if len(values) != len(set(values)):
        raise argparse.ArgumentTypeError("duplicate values are not allowed")
    return values


def save(path: Path, report: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n")
    temporary.replace(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--presets", nargs="+", choices=PRESETS, default=list(PRESETS))
    parser.add_argument("--counts", type=lambda value: parse_csv(value, positive=True),
                        default=list(STANDARD_DURATIONS))
    parser.add_argument("--seeds", type=lambda value: parse_csv(value), default=[123])
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--duration", type=float,
                        help="override the standard count-dependent physical duration")
    parser.add_argument("--shaft-speed", type=float, default=.15)
    parser.add_argument("--numeric-backend", choices=("auto", "python", "numba"),
                        default="auto")
    parser.add_argument("--pair-search", choices=("grid", "sweep", "all"),
                        default="sweep")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.repeats < 1 or args.duration is not None and args.duration <= 0:
        parser.error("repeats and duration must be positive")
    if args.shaft_speed <= 0:
        parser.error("shaft speed must be positive")
    unknown = [count for count in args.counts if count not in STANDARD_DURATIONS]
    if unknown and args.duration is None:
        parser.error("nonstandard counts require --duration")

    report = {
        "schema_version": 1,
        "complete": False,
        "python_executable": sys.executable,
        "python": platform.python_version(),
        "settings": {"presets": args.presets, "counts": args.counts,
                     "seeds": args.seeds, "repeats": args.repeats,
                     "duration_override": args.duration,
                     "standard_durations": STANDARD_DURATIONS,
                     "shaft_speed": args.shaft_speed,
                     "numeric_backend": args.numeric_backend,
                     "pair_search": args.pair_search},
        "runs": [],
    }
    save(args.output, report)
    for preset in args.presets:
        for count in args.counts:
            for seed in args.seeds:
                for repeat in range(args.repeats):
                    duration = args.duration or STANDARD_DURATIONS[count]
                    command = [sys.executable, "-m",
                               "microthermo.measurements.performance_benchmark",
                               "--preset", preset, "--particles", str(count),
                               "--duration", str(duration), "--seed", str(seed),
                               "--shaft-speed", str(args.shaft_speed),
                               "--pair-search", args.pair_search,
                               "--numeric-backend", args.numeric_backend]
                    start = time.perf_counter()
                    result = subprocess.run(command, capture_output=True, text=True)
                    elapsed = time.perf_counter() - start
                    if result.returncode:
                        report["runs"].append({"preset": preset, "particles": count,
                                               "seed": seed, "repeat": repeat,
                                               "error": result.stderr.strip(),
                                               "process_seconds": elapsed})
                        save(args.output, report)
                        print(result.stderr, file=sys.stderr)
                        return result.returncode
                    row = json.loads(result.stdout)
                    row["repeat"] = repeat
                    row["process_seconds"] = elapsed
                    report["runs"].append(row)
                    save(args.output, report)
                    print(f"{preset} n={count} seed={seed} repeat={repeat}: "
                          f"{row['physical_per_wall_second']:.3g} physical s/wall s",
                          flush=True)
    report["complete"] = True
    save(args.output, report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
