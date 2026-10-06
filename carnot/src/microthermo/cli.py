from __future__ import annotations

import argparse
import json
import math
import platform
from pathlib import Path
import sys
import time
import tomllib
import numpy as np

from . import __version__
from .api import load_preset
from .config import RunConfig
from .experiments import preset_names
from .io.exports import export_run
from .io.replay import precalculate_replay
from .measurements.speed_study import run_speed_study
from .validation import results_as_dict, run_validation


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="microthermo")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="run a reproducible headless experiment")
    run.add_argument("--config", help="TOML experiment configuration")
    run.add_argument("--preset", choices=preset_names())
    run.add_argument("--seed", type=int)
    run.add_argument("--duration", type=float)
    run.add_argument("--cycles", type=int)
    run.add_argument("--particles", type=int)
    run.add_argument("--temperature", type=float)
    run.add_argument("--carnot-radius-scale", type=float,
                     help="multiply Carnot particle radius; default 1.0")
    run.add_argument("--sample-interval", type=float)
    run.add_argument("--max-horizon", type=float)
    run.add_argument("--pair-search", choices=("grid", "sweep", "all"))
    run.add_argument("--engine", choices=("reference", "edmd"),
                     help="edmd: compiled event-driven kernel (default for osmosis presets)")
    _osmosis_arguments(run)
    run.add_argument("--numeric-backend", choices=("auto", "python", "numba"))
    run.add_argument("--wall-search", choices=("bounded", "all"))
    run.add_argument("--wall-kernel", choices=("auto", "python"))
    run.add_argument("--penetration-kernel", choices=("auto", "python"))
    run.add_argument("--pair-kernel", choices=("auto", "scalar"))
    run.add_argument("--cam-kernel", choices=("auto", "python"))
    run.add_argument("--output", default="run")
    run.add_argument("--headless", action="store_true", help="accepted for command compatibility")
    run.add_argument("--reversed", action="store_true")
    run.add_argument("--shaft-mode", choices=("controlled", "free"))
    run.add_argument("--shaft-speed", type=float,
                     help="Carnot shaft angular speed in rad/s (default: 0.15)")
    run.add_argument("--cold-jacket", action="store_true", default=None,
                     help="add cold-sector thermal top and bottom walls")
    run.add_argument("--hot-jacket", action="store_true", default=None,
                     help="add hot-sector thermal top and bottom walls")
    run.add_argument("--cam-fractions", type=float, nargs=4,
                     metavar=("HOT", "EXPANSION", "COLD", "COMPRESSION"),
                     help="fractions of a revolution assigned to the four cam sectors")
    run.add_argument("--initial-temperature", type=float,
                     help="Carnot: starting gas temperature (default: --temperature)")
    run.add_argument("--temperature-ratio", type=float,
                     help="Carnot: T_hot/T_cold with T_hot = 1.5*temperature (default 2.0)")
    run.add_argument("--transient-cycles", type=int)
    run.add_argument("--efficiency-min-cycles", type=int)
    study = sub.add_parser("speed-study", help="compare controlled Carnot shaft speeds across seeds")
    study.add_argument("--preset", choices=("carnot_discs","carnot_triangles"),
                       default="carnot_discs")
    study.add_argument("--speeds", type=float, nargs="+", default=[.075,.15,.3])
    study.add_argument("--seeds", type=int, nargs="+", default=[123,124,125])
    study.add_argument("--cycles", type=int, default=10)
    study.add_argument("--particles", type=int, default=8)
    study.add_argument("--carnot-radius-scale", type=float, default=1.0,
                       help="multiply Carnot particle radius; default 1.0")
    study.add_argument("--max-horizon", type=float, default=.05)
    study.add_argument("--transient-cycles", type=int, default=2)
    study.add_argument("--efficiency-min-cycles", type=int, default=8)
    study.add_argument("--cold-jacket", action="store_true",
                       help="add cold-sector thermal top and bottom walls")
    study.add_argument("--hot-jacket", action="store_true",
                       help="add hot-sector thermal top and bottom walls")
    study.add_argument("--cam-fractions", type=float, nargs=4,
                       metavar=("HOT", "EXPANSION", "COLD", "COMPRESSION"),
                       default=(.25,.25,.25,.25),
                       help="fractions of a revolution assigned to the four cam sectors")
    study.add_argument("--output", default="runs/carnot_speed_study.json")
    replay = sub.add_parser("precalculate", help="write chunked replay frames offline")
    replay.add_argument("--preset", choices=preset_names(), default="carnot_discs")
    replay.add_argument("--seed", type=int, default=123)
    replay.add_argument("--engine", choices=("reference", "edmd"))
    _osmosis_arguments(replay)
    replay.add_argument("--particles", type=int, default=32)
    length = replay.add_mutually_exclusive_group(required=True)
    length.add_argument("--duration", type=float)
    length.add_argument("--cycles", type=int, help="Carnot: whole shaft revolutions")
    replay.add_argument("--shaft-speed", type=float, default=.15)
    replay.add_argument("--cold-jacket", action="store_true")
    replay.add_argument("--hot-jacket", action="store_true")
    replay.add_argument("--cam-fractions", type=float, nargs=4,
                        metavar=("HOT", "EXPANSION", "COLD", "COMPRESSION"),
                        default=(.25,.25,.25,.25))
    replay.add_argument("--initial-temperature", type=float,
                     help="Carnot: starting gas temperature (default: --temperature)")
    replay.add_argument("--temperature-ratio", type=float,
                     help="Carnot: T_hot/T_cold with T_hot = 1.5*temperature (default 2.0)")
    replay.add_argument("--pair-search", choices=("grid", "sweep", "all"), default="sweep",
                        help="broad phase; sweep is about 2x faster for 500 triangles")
    replay.add_argument("--fps", type=float, default=60.)
    replay.add_argument("--chunk-frames", type=int, default=256)
    replay.add_argument("--output", required=True)
    replay_view = sub.add_parser("replay", help="play a precalculated archive")
    replay_view.add_argument("path")
    val = sub.add_parser("validate", help="run physics validation gates")
    val.add_argument("--suite", default="scientific", choices=("scientific", "quick"))
    bench = sub.add_parser("benchmark", help="measure reference engine throughput")
    bench.add_argument("--suite", default="standard")
    gui = sub.add_parser("gui", help="launch the optional Qt laboratory")
    gui.add_argument("--preset", choices=preset_names(), default="gas_box")
    gui.add_argument("--transient-cycles", type=int, default=2)
    gui.add_argument("--efficiency-min-cycles", type=int, default=8)
    gui.add_argument("--shaft-speed", type=float, default=.15)
    gui.add_argument("--shaft-mode", choices=("controlled", "free"), default="controlled")
    gui.add_argument("--seed", type=int, default=123)
    gui.add_argument("--engine", choices=("reference", "edmd"))
    _osmosis_arguments(gui)
    gui.add_argument("--particles", type=int,
                     help="default: 96 for carnot_triangles, otherwise 48")
    gui.add_argument("--cold-jacket", action="store_true")
    gui.add_argument("--hot-jacket", action="store_true")
    gui.add_argument("--initial-temperature", type=float,
                     help="Carnot: starting gas temperature (default: --temperature)")
    gui.add_argument("--temperature-ratio", type=float, default=2.0,
                     help="Carnot: T_hot/T_cold with T_hot = 1.5*temperature (default 2.0)")
    return p


_OSMOSIS_KEYS = {"hosts": 8, "binding_energy": None, "host_mass": 25.0,
                 "host_area_fraction": 0.2, "mouth_width": None, "rough_fraction": 1.0,
                 "discs_start": "right", "temperature_schedule": ()}


def _temperature_steps(text: str) -> tuple[tuple[float, float], ...]:
    """Parse "100:2,300:0.75" into ((100., 2.), (300., .75))."""
    try:
        return tuple((float(t), float(v)) for t, v in
                     (item.split(":") for item in text.split(",") if item.strip()))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("use time:temperature pairs, e.g. 100:2,300:0.75") from exc


def _osmosis_arguments(sub) -> None:
    group = sub.add_argument_group("osmosis presets")
    group.add_argument("--hosts", type=int, help="number of ring hosts (default 8)")
    group.add_argument("--binding-energy", type=float,
                       help="well depth epsilon (default 1.5; 0 for osmosis_hard)")
    group.add_argument("--host-mass", type=float, help="host mass (default 25)")
    group.add_argument("--host-area-fraction", type=float,
                       help="fraction of the left chamber covered by hosts (default 0.2)")
    group.add_argument("--mouth-width", type=float,
                       help="clear mouth chord (default three disc diameters)")
    group.add_argument("--rough-fraction", type=float,
                       help="probability that a host contact is rough (default 1)")
    group.add_argument("--discs-start", choices=("right", "mixed"),
                       help="where the discs start (default right of the membrane)")
    group.add_argument("--temperature-steps", dest="temperature_schedule",
                       type=_temperature_steps,
                       help="wall temperature changes as time:T pairs, e.g. 300:2,600:0.75")


def _osmosis_values(args) -> dict:
    return {key: getattr(args, key) for key in _OSMOSIS_KEYS
            if getattr(args, key) is not None}


def _run(args) -> int:
    values={"preset":"gas_box","seed":123,"duration":10.,"particles":None,
            "temperature":1.,"sample_interval":.05,"max_horizon":.05,
            "cycles":None,"reversed_cycle":False,"shaft_mode":"controlled",
            "shaft_speed":.15,
            "carnot_radius_scale":1.0,
            "transient_cycles":2,"efficiency_min_cycles":8,
            "pair_search":"grid","numeric_backend":"auto",
            "wall_search":"bounded","wall_kernel":"auto",
            "penetration_kernel":"auto","pair_kernel":"auto",
            "cam_kernel":"auto","initial_temperature":None,
            "temperature_ratio":2.0,"engine":None,**_OSMOSIS_KEYS}
    if args.config:
        with open(args.config,"rb") as f:
            loaded=tomllib.load(f)
        loaded=loaded.get("experiment",loaded)
        unknown=set(loaded)-set(values)
        if unknown: raise SystemExit(f"unknown configuration keys: {', '.join(sorted(unknown))}")
        values.update(loaded)
    for key in ("preset","seed","duration","particles","temperature","carnot_radius_scale","sample_interval",
                "max_horizon","cycles","shaft_mode","shaft_speed","transient_cycles",
                "efficiency_min_cycles","pair_search","numeric_backend",
                "wall_search","wall_kernel","penetration_kernel","pair_kernel",
                "cam_kernel","cold_jacket","hot_jacket","cam_fractions",
                "initial_temperature","temperature_ratio","engine",*_OSMOSIS_KEYS):
        value=getattr(args,key)
        if value is not None: values[key]=tuple(value) if key=="cam_fractions" else value
    if args.reversed: values["reversed_cycle"]=True
    duration = values["duration"]
    if values["cycles"] is not None:
        if not values["preset"].startswith("carnot"):
            raise SystemExit("--cycles applies to Carnot presets")
        if values["cycles"] <= 0:
            raise SystemExit("--cycles must be positive")
        if not math.isfinite(values["shaft_speed"]) or values["shaft_speed"] <= 0:
            raise SystemExit("shaft speed must be finite and positive")
        duration = values["cycles"]*2*math.pi/values["shaft_speed"]
    values["duration"]=duration
    cfg = RunConfig(**values)
    sim = load_preset(cfg)
    next_sample = min(cfg.sample_interval, duration)
    while sim.time < duration-1e-12:
        sim.samples.append(sim.advance_to(min(next_sample, duration)))
        next_sample = min(duration, next_sample+cfg.sample_interval)
    output = export_run(args.output, cfg, sim)
    print(json.dumps({"output": str(output), "events": sim.event_count,
                      "energy_residual": sim.snapshot().energy_residual}, indent=2))
    return 0


def _validate(args) -> int:
    results = run_validation(args.suite)
    print(json.dumps(results_as_dict(results), indent=2))
    return 0 if all(r.passed for r in results) else 2


def _benchmark(args) -> int:
    rows=[]
    # This measures the auditable all-pairs reference scheduler. Large-scene
    # benchmark matrices belong to the future compiled broad-phase backend.
    for n in (20, 50):
        physical=.1
        cfg=RunConfig("gas_box", 123, physical, n, 1., physical, .05)
        sim=load_preset(cfg); start=time.perf_counter(); sim.advance_to(physical); elapsed=time.perf_counter()-start
        rows.append({"particles":n,"physical_seconds":physical,"wall_seconds":elapsed,
                     "events":sim.event_count,"events_per_second":sim.event_count/elapsed})
    print(json.dumps({"platform":platform.platform(),"python":platform.python_version(),
                      "numpy":np.__version__,"results":rows},indent=2))
    return 0


def _speed_study(args) -> int:
    try:
        result=run_speed_study(
            preset=args.preset,speeds=tuple(args.speeds),seeds=tuple(args.seeds),
            cycles=args.cycles,particles=args.particles,
            max_horizon=args.max_horizon,
            transient_cycles=args.transient_cycles,
            efficiency_min_cycles=args.efficiency_min_cycles,
            cold_jacket=args.cold_jacket,hot_jacket=args.hot_jacket,
            cam_fractions=tuple(args.cam_fractions),
            carnot_radius_scale=args.carnot_radius_scale)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    output=Path(args.output)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"output":str(output),"speeds":[
        {"shaft_speed":group["shaft_speed"],
         "ready_seed_count":group["ready_seed_count"],
         "mean_measured_efficiency":group["mean_measured_efficiency"]}
        for group in result["speeds"]]},indent=2))
    return 0


def main(argv=None) -> int:
    args=parser().parse_args(argv)
    if args.command == "run": return _run(args)
    if args.command == "validate": return _validate(args)
    if args.command == "benchmark": return _benchmark(args)
    if args.command == "speed-study": return _speed_study(args)
    if args.command == "precalculate":
        duration=args.duration
        if args.cycles is not None:
            if args.cycles <= 0 or not args.preset.startswith("carnot_"):
                raise SystemExit("--cycles needs a Carnot preset and a positive count")
            duration=args.cycles*2*math.pi/args.shaft_speed
        cfg=RunConfig(preset=args.preset,seed=args.seed,particles=args.particles,
                      duration=duration,shaft_speed=args.shaft_speed,
                      cold_jacket=args.cold_jacket,hot_jacket=args.hot_jacket,
                      cam_fractions=tuple(args.cam_fractions),
                      initial_temperature=args.initial_temperature,
                      temperature_ratio=(2.0 if args.temperature_ratio is None
                                         else args.temperature_ratio),
                      pair_search=args.pair_search,engine=args.engine,
                      **_osmosis_values(args))
        started=time.perf_counter()
        def progress(physical,duration,events):
            elapsed=time.perf_counter()-started
            remaining=elapsed*(duration-physical)/physical if physical>0 else math.nan
            print(f"t={physical:.1f}/{duration:.1f} s ({100*physical/duration:.1f}%) "
                  f"events={events} elapsed={elapsed/60:.1f} min "
                  f"remaining≈{remaining/60:.0f} min",file=sys.stderr,flush=True)
        output=precalculate_replay(cfg,args.output,fps=args.fps,
                                   chunk_frames=args.chunk_frames,progress=progress)
        print(json.dumps({"output":str(output)},indent=2))
        return 0
    if args.command == "replay":
        try:
            from .ui.replay_window import launch_replay
        except ImportError as exc:
            print("GUI dependencies are missing; install with: pip install -e '.[gui]'", file=sys.stderr)
            return 3
        return launch_replay(args.path)
    if args.command == "gui":
        if args.transient_cycles < 0 or args.efficiency_min_cycles < 4:
            raise SystemExit("require nonnegative transients and at least four efficiency cycles")
        if not math.isfinite(args.shaft_speed) or args.shaft_speed <= 0:
            raise SystemExit("shaft speed must be finite and positive")
        try:
            from .ui.main_window import launch
        except ImportError as exc:
            print("GUI dependencies are missing; install with: pip install -e '.[gui]'", file=sys.stderr)
            return 3
        config=RunConfig(preset=args.preset,seed=args.seed,
                         particles=(args.particles if args.particles is not None else
                                    96 if args.preset=="carnot_triangles" else
                                    200 if args.preset.startswith("osmosis") else 48),
                         max_horizon=.02,shaft_speed=args.shaft_speed,
                         pair_search="sweep",
                         shaft_mode=args.shaft_mode,cold_jacket=args.cold_jacket,
                         hot_jacket=args.hot_jacket,
                         initial_temperature=args.initial_temperature,
                         temperature_ratio=args.temperature_ratio,
                         transient_cycles=args.transient_cycles,
                         efficiency_min_cycles=args.efficiency_min_cycles,
                         engine=args.engine,**_osmosis_values(args))
        try:
            config.validate()
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
        return launch(args.preset,args.transient_cycles,args.efficiency_min_cycles,
                      args.shaft_speed,args.shaft_mode,config=config)
    return 1
