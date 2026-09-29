# Microscopic Thermodynamics Laboratory

This repository contains a deterministic, event-driven reference model for the
two experiments in `microscopic_thermodynamics_implementation.md`: a
particle-driven Carnot apparatus and a labyrinth association laboratory.

The model uses normalized two-dimensional units (`k_B = 1`).  Hard contacts
are resolved from impulses, thermal walls refresh the contact-normal mode, and
heat and work are recorded event by event.  The reference backend emphasizes
auditability.  Optional Numba and Qt dependencies are separated from the
headless package.

## Install and run

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pip install -e .
```

For the exact headless dependency tested here, install
`-r requirements-lock.txt` before the editable package.

Install `.[gui]` to use `python -m microthermo gui`.  Headless commands never
import Qt.  Output directories contain the resolved configuration, summary,
sampled CSV, and event JSONL.

The default scenes are intentionally modest.  The rotating polygon CCD and
compound-host machinery are strict reference implementations suitable for
validation and extension; large production scenes should use the optional
compiled kernels after profiling.

`carnot_triangles` remains experimental while broader collision validation,
mixed-policy contact clusters, and the mechanism/UI acceptance gates are
completed. Numerical failures pause the GUI and write a diagnostic checkpoint
in a temporary directory. See
`docs/CARNOT_GAP_CLOSURE_PLAN.md` for the remaining acceptance gates.

## Available configurations

Run these commands from the repository root after activating the virtual
environment. Each command uses a separate output directory so results are not
mixed between experiments.

### Gas and equilibration demonstrations

Thermalized spinless disc gas:

```bash
python -m microthermo run \
  --preset gas_box \
  --duration 20 \
  --particles 32 \
  --temperature 1.0 \
  --seed 123 \
  --output runs/gas_box
```

Rotating triangles relaxing toward translational and rotational
equipartition:

```bash
python -m microthermo run \
  --preset triangle_equipartition \
  --duration 20 \
  --particles 32 \
  --temperature 1.0 \
  --seed 123 \
  --output runs/triangle_equipartition
```

### Carnot apparatus

Controlled-speed disc engine:

```bash
python -m microthermo run \
  --preset carnot_discs \
  --shaft-mode controlled \
  --cycles 1 \
  --particles 32 \
  --seed 123 \
  --output runs/carnot_discs_controlled
```

Free-running disc engine:

```bash
python -m microthermo run \
  --preset carnot_discs \
  --shaft-mode free \
  --cycles 1 \
  --particles 32 \
  --max-horizon 0.005 \
  --seed 123 \
  --output runs/carnot_discs_free
```

Controlled-speed triangle engine:

```bash
python -m microthermo run \
  --preset carnot_triangles \
  --shaft-mode controlled \
  --cycles 1 \
  --particles 32 \
  --max-horizon 0.01 \
  --seed 123 \
  --output runs/carnot_triangles_controlled
```

Free-running triangle engine:

```bash
python -m microthermo run \
  --preset carnot_triangles \
  --shaft-mode free \
  --cycles 1 \
  --particles 32 \
  --max-horizon 0.005 \
  --seed 123 \
  --output runs/carnot_triangles_free
```

Reversed controlled disc cycle, operating as a refrigerator or heat pump:

```bash
python -m microthermo run \
  --preset carnot_discs \
  --shaft-mode controlled \
  --reversed \
  --cycles 1 \
  --particles 32 \
  --seed 123 \
  --output runs/carnot_discs_reversed
```

The `--reversed` option can also be combined with `carnot_triangles` and with
`--shaft-mode free`. The `--cycles` option applies only to Carnot presets; use
`--duration` for the other experiments. `--shaft-speed` sets the physical
angular speed in rad/s (default `0.15`). In controlled mode, `--cycles` uses
that speed to stop at the requested number of full shaft turns. In free mode,
it sets only a nominal duration; the shaft may not complete that many turns.

To request a measured net efficiency, run at least the transient count plus
the minimum measurement count in complete cycles. The defaults skip two cycles
and require eight more:

```bash
python -m microthermo run \
  --preset carnot_discs --shaft-mode controlled \
  --cycles 10 --particles 32 --seed 123 \
  --transient-cycles 2 --efficiency-min-cycles 8 \
  --output runs/carnot_discs_efficiency
```

The same two gate options apply to `microthermo gui`. `summary.json` records
the gate status and reports `efficiency` only when the matched cycles have
safely positive hot heat in every contiguous uncertainty block. The estimate
is signed `(load output − motor work) / hot heat`; its four-block jackknife
standard error describes between-block variation, not a confidence interval.
The cycle count and total stored-energy change accompany the estimate. Small
or fast runs can legitimately report no efficiency.

To collect a reproducible speed comparison across independent seeds:

```bash
python -m microthermo speed-study \
  --preset carnot_discs --speeds 0.075 0.15 0.3 \
  --seeds 123 124 125 --cycles 10 --particles 8 \
  --output runs/carnot_speed_study.json
```

The JSON contains completed-cycle accounting, each seed's efficiency gate,
branch-by-branch thermal contact counts, positive and negative heat, boundary
temperatures and gas energy, descriptive early/late gas-energy and hot-heat
drift screens, and the analytical ideal reference. Thermal event records also
include incoming and outgoing contact-normal mode energies. Their difference is
reconciled with the existing heat ledger; these contact energies are distinct
from the whole-gas translational temperature. A speed's mean measured
efficiency and seed standard error appear only when at least two seeds pass
the hot-heat gate.
The drift screens require at least 16 cycles after transients. Each compares
four contiguous cycle-batch means from the early and late halves against a
10% magnitude tolerance, with a two-standard-error descriptive margin.
`shift_bounded` means only that this measured early/late shift fits that
tolerance; batches may remain correlated, so it does not certify stationarity.
These are descriptive statistics, not a convergence verdict. A slow speed
requires more physical simulation time and can take substantially longer.

For the minimum cycle-batch stationarity screen with default two transient
cycles, request at least 18 cycles:

```bash
python -m microthermo speed-study \
  --preset carnot_discs --speeds 0.15 \
  --seeds 123 124 --cycles 18 --particles 16 \
  --output runs/carnot_stationarity.json
```

### Chemical-potential and permeability demonstrations

Hard-wall labyrinth with geometric association:

```bash
python -m microthermo run \
  --preset labyrinth \
  --duration 10 \
  --particles 32 \
  --temperature 1.0 \
  --seed 123 \
  --output runs/labyrinth_hard
```

Labyrinth with a reversible interior binding-energy step:

```bash
python -m microthermo run \
  --preset labyrinth_energetic \
  --duration 10 \
  --particles 32 \
  --temperature 1.0 \
  --seed 123 \
  --output runs/labyrinth_energetic
```

Repeated-gap membrane that passes small discs and rejects large mobile hosts:

```bash
python -m microthermo run \
  --preset selective_membrane \
  --duration 10 \
  --particles 32 \
  --temperature 1.0 \
  --seed 123 \
  --output runs/selective_membrane
```

### Configuration file

Run the included free-shaft Carnot configuration:

```bash
python -m microthermo run \
  --config examples/carnot_free.toml \
  --output runs/carnot_from_toml
```

Command-line options override values from the TOML file. For example:

```bash
python -m microthermo run \
  --config examples/carnot_free.toml \
  --duration 2 \
  --particles 32 \
  --output runs/carnot_toml_short
```

### Validation and benchmark commands

```bash
python -m microthermo validate --suite scientific
python -m microthermo validate --suite quick
python -m unittest discover -s tests -v
python -m microthermo benchmark --suite standard
```

### Optional graphical interface

Install the GUI dependencies and launch any preset:

```bash
python -m pip install -e '.[gui]'
python -m microthermo gui --preset gas_box
python -m microthermo gui --preset triangle_equipartition
python -m microthermo gui --preset carnot_discs
python -m microthermo gui --preset carnot_triangles
python -m microthermo gui --preset labyrinth
python -m microthermo gui --preset labyrinth_energetic
python -m microthermo gui --preset selective_membrane
```

In the instrument pane, open **Numerical diagnostics** to see maximum
penetration, CCD refinements and failures, contact clusters and their maximum
solver residual, event throughput, and achieved playback rate. Throughput is
events per wall-clock second; playback rate is physical simulation seconds per
wall-clock second. Both rates show **Paused** outside continuous playback.

Use smaller particle counts for quick reference-backend runs. The scheduler
prioritizes auditable collision chronology over large-scene speed.

For larger Carnot gases, particle radius scales as `0.025*sqrt(48/particles)`.
Triangles occupy about 3.12% and discs about 7.54% of the minimum cylinder
area at every particle count. The default swept grid searches for possible
body pairs; the optional Numba package compiles its numeric reach filter and
polygon pair collision search.

```bash
.venv/bin/python -m pip install -r requirements-accel-lock.txt
.venv/bin/python -m pip install -e '.[accel]'
.venv/bin/python -m microthermo run --preset carnot_triangles \
  --particles 128 --duration 0.1 --pair-search grid \
  --numeric-backend numba --wall-search bounded \
  --wall-kernel auto --penetration-kernel auto \
  --output runs/triangles_128
.venv/bin/python -m microthermo.measurements.performance_benchmark \
  --preset carnot_triangles --particles 64 --duration 0.1 \
  --pair-search grid --numeric-backend numba
```

For the unpruned reference, use `--pair-search all --numeric-backend python`
with `--wall-search all --wall-kernel python --penetration-kernel python`.
Numba is optional; `--numeric-backend auto` uses it when installed and falls
back to Python otherwise. Short 32–64 particle triangle runs are about 8×
faster with the compiled geometry kernel on the measured machine. The full
compiled pair collision loop, bounded wall checks, and wall kernels improve
short 200/500 triangle probes further. They remain far below interactive
playback, and long-run physics still needs validation. See the
[performance plan](docs/CARNOT_GAP_CLOSURE_PLAN.md).

Example TOML:

```toml
[experiment]
preset = "carnot_triangles"
seed = 123
duration = 20.0
particles = 64
temperature = 1.0
sample_interval = 0.05
max_horizon = 0.01
shaft_mode = "free"
```

## Implemented scope

- Exact ballistic flight, analytic disc TOI, finite segment walls with polygon
  contact geometry, rotating convex polygon conservative advancement, elastic
  rotational impulses, small stationary elastic contact clusters, and strict
  failure on unresolved collision queries.
- Specular and flux-distributed thermal contacts with event-level heat, work,
  support impulse, and compensated energy ledgers.
- C2 four-branch Carnot cam, shared physical selector, controlled and
  finite-inertia free shaft modes with reflected piston mass, load accounting,
  reversed operation, disc and triangle gases. Carnot snapshots include named
  apparatus component paths and energy state; the Qt preview draws the named
  apparatus with a fitted, aspect-preserving camera.
- Hard and energetic nested labyrinths with physical apertures, trajectory
  portal events, occupancy labels, and reversible square-well impulses.
- A repeated-gap selective membrane with mobile large host particles and small
  permeating discs.
- Deterministic checkpoints, JSON/CSV event exports, scientific validation,
  benchmarks, and an optional Qt view that receives read-only snapshots from
  a simulation worker thread. The Qt view includes a fitted apparatus,
  pressure–area and temperature plots, and a scrollable energy ledger.
  Completed Carnot cycles include exact-boundary gas/apparatus ledger summaries
  in `cycles.json`. For controlled-speed cycles, the live panel and export also
  compare event-level gas work with the pressure–area loop integral. The
  pressure instrument uses a fixed 0.25 physical-time trailing impulse window
  (zero-padded at startup); the reported smoothing bound covers the difference
  caused by averaging and cycle-boundary clipping. Free-shaft cycles have no
  pressure–area work comparison yet. Dashed curves on the pressure–area plot
  are an analytical point-gas Carnot reference: `P = N T / A` on isotherms and
  `T A^(2/f) = constant` on adiabats, with `f = 2` for smooth discs and `f = 3`
  for rotating triangles. The label gives the theoretical `ηC = 1 − Tc/Th`;
  it is not an efficiency measured from the run.

The current labyrinth walls are fixed finite segments.  Mobile compound hosts,
per-host BVHs, general mixed-policy contact clusters, spring-supported
porous partitions, finite particle reservoirs, further compiled Numba kernels, and the
full polished plotting/recording interface remain extension work.  The engine
stops on ambiguous event chronology; it does not claim validated dense or
jammed multi-contact dynamics.
