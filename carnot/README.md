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

The desktop toolbar has play/pause, collision, duration, and Carnot branch
steps; seed reset; preset selection; and separate target playback and physical
shaft-speed controls. A shaft-speed edit is recorded as a motor intervention.
Use the mouse wheel to zoom, drag empty scene space to pan, click a particle
to inspect it, and select **Fit scene** to restore the view.

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

When `--particles` is omitted, `carnot_triangles` now starts with 96
triangles; an explicit count still takes precedence. The GUI uses the same
96-triangle default for this preset.

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

For the larger thermodynamic evidence matrices in the gap plan:

```bash
.venv/bin/python -m microthermo speed-study \
  --preset carnot_discs --speeds 0.3 0.15 0.075 \
  --seeds 123 124 125 126 127 --cycles 18 --particles 16 \
  --output runs/carnot_thermo_16disc_5seed_3speed_18cycles.json

.venv/bin/python -m microthermo speed-study \
  --preset carnot_discs --speeds 0.075 \
  --seeds 123 124 125 --cycles 18 --particles 32 \
  --output runs/carnot_thermo_32disc_3seed_slow_18cycles.json

.venv/bin/python -m microthermo speed-study \
  --preset carnot_triangles --speeds 0.15 0.075 \
  --seeds 123 124 125 --cycles 18 --particles 32 \
  --output runs/carnot_thermo_32triangle_3seed_2speed_18cycles.json
```

The speed study reports descriptive hot-heat and gas-energy drift screens.
An efficiency estimate requires positive hot heat in each contiguous block;
it does not by itself establish stationarity or convergence to the ideal
Carnot value. The measured results and limitations of these matrices are in
[`docs/CARNOT_THERMODYNAMIC_EVIDENCE_2026-09-30.md`](docs/CARNOT_THERMODYNAMIC_EVIDENCE_2026-09-30.md).

An optional cold jacket makes the stationary top and bottom cylinder walls
thermal at the cold-reservoir temperature during the cold cam sector. They
remain specular in the other sectors. The default apparatus has no jacket.
This is an experimental change to reservoir coupling, not a temperature clamp.
For a matched controlled-shaft comparison, add `--cold-jacket` to a `run` or
`speed-study` command (or set `cold_jacket = true` in a run TOML file):

```bash
.venv/bin/python -m microthermo speed-study \
  --preset carnot_triangles --speeds 0.075 \
  --seeds 123 124 125 --cycles 18 --particles 32 --cold-jacket \
  --output runs/carnot_cold_jacket_32triangle.json
```

The cold-jacket comparison and its limits are recorded in the thermodynamic
evidence document linked above. The efficiency and drift gates are unchanged.

To increase hot-sector contact area while retaining the cold jacket, add
`--hot-jacket` alongside `--cold-jacket` (or set `hot_jacket = true` in a run
TOML file). The top and bottom walls are thermal at the hot temperature only
in the hot sector, thermal at the cold temperature only in the cold sector,
and specular in the two adiabatic sectors. The selector wall and piston cam
are unchanged. This is a separate apparatus protocol; the same efficiency
and drift gates apply.

To change cam sector timing, pass four positive fractions in hot, adiabatic
expansion, cold, adiabatic compression order. They must sum to one. For
example, this protocol keeps hot and cold exposure at one quarter-cycle each,
shortens expansion to 15%, and lengthens compression to 35%:

```bash
.venv/bin/python -m microthermo speed-study \
  --preset carnot_discs --speeds 0.075 \
  --seeds 123 124 125 --cycles 18 --particles 16 --cold-jacket \
  --cam-fractions 0.25 0.15 0.25 0.35 \
  --output runs/carnot_long_compression.json
```

The same `--cam-fractions` option works for `run`; a run TOML file may set
`cam_fractions = [0.25, 0.15, 0.25, 0.35]`. The default is four equal
sectors. Timing changes preserve the cam's area endpoints but change piston
speed within the affected sectors.

To repeat the longer, three-seed slow triangle jacket study:

```bash
.venv/bin/python -m microthermo speed-study \
  --preset carnot_triangles --speeds 0.075 \
  --seeds 123 124 125 --cycles 36 --particles 32 --cold-jacket \
  --output runs/carnot_cold_jacket_32triangle_3seed_36cycles.json
```

This took about 38 minutes on an i7-11850H. It improved positive-hot-heat
readiness but did not establish stationarity.

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

The Carnot canvas labels the default controlled-speed shaft as a motor, since
its angular speed is imposed. The circular housing marks its phase; it is not
an independently moving flywheel. In free-shaft runs the same coordinate is
shown as a flywheel with finite inertia. The stylized rotating spring coil is
hidden; an active spring is named with its stored energy, which is also shown
in the instruments. The cam and thermal
selector depict actual piston and contact rules without separate dynamic
coordinates. The load label displays accumulated work, not a moving storage
body.

To compare the two shaft modes in the desktop view:

```bash
.venv/bin/python -m microthermo gui --preset carnot_triangles --shaft-mode controlled
.venv/bin/python -m microthermo gui --preset carnot_triangles --shaft-mode free
```

Use smaller particle counts for quick reference-backend runs. The scheduler
prioritizes auditable collision chronology over large-scene speed.

For larger Carnot gases, particle radius scales as `0.025*sqrt(48/particles)`.
Triangles occupy about 3.12% and discs about 7.54% of the minimum cylinder
area at every particle count. The default swept grid searches for possible
body pairs; the optional Numba package compiles its numeric reach filter and
polygon pair collision search. For larger runs, `--pair-search sweep` selects
a compiled swept-box search when Numba is installed.
The default `--pair-kernel auto` batches triangle pair collision queries;
it also batches disc pair queries. `--pair-kernel scalar` keeps the scalar
query path for comparison. Numba also batches conservative wall reach checks
and controlled-speed disc contacts with the cam piston; `--cam-kernel python`
selects the cam reference path. Controlled-speed triangle contacts with the
cam piston are also batched. Free-shaft piston motion remains on its
general trajectory solver. Fixed-wall disc contacts also use a compiled batch;
`--wall-kernel python` selects their scalar reference path. Disc-only overlap
checks also use a compiled batch in infinite-wall worlds;
`--penetration-kernel python` selects the reference path.
The batched disc scheduler decodes contacts only near the next event time,
while retaining all failure checks and simultaneous contacts.
Triangle-only worlds also batch current-state pair overlap checks; the same
`--penetration-kernel python` option keeps the scalar diagnostic path.
Checkpoint copies share the immutable cam geometry and copy the mutable shaft
and compensated ledger counters. Energy checks skip region accounting in
worlds without portals.

```bash
.venv/bin/python -m pip install -r requirements-accel-lock.txt
.venv/bin/python -m pip install -e '.[accel]'
.venv/bin/python -m microthermo run --preset carnot_triangles \
  --particles 128 --duration 0.1 --pair-search sweep \
  --numeric-backend numba --pair-kernel auto --wall-search bounded \
  --wall-kernel auto --penetration-kernel auto \
  --output runs/triangles_128
.venv/bin/python -m microthermo.measurements.performance_benchmark \
  --preset carnot_triangles --particles 200 --duration 0.02 \
  --pair-search sweep --numeric-backend numba
```

For the unpruned reference, use `--pair-search all --numeric-backend python`
with `--pair-kernel scalar --wall-search all --wall-kernel python`
and `--penetration-kernel python --cam-kernel python`.
Numba is optional; `--numeric-backend auto` uses it when installed and falls
back to Python otherwise. Short 32–64 particle triangle runs are about 8×
faster with the compiled geometry kernel on the measured machine. The full
compiled pair collision loop, bounded wall checks, and wall kernels improve
short 200/500 triangle probes further. The optional swept-box search improved
short 200 and 500 triangle probes over the grid on the measured machine. They
remain far below interactive playback, and long-run physics still needs
validation. A compressed offline replay archive and desktop playback are
available; long-run display acceptance remains open. See the
[performance plan](docs/CARNOT_GAP_CLOSURE_PLAN.md).

To precalculate a versioned, chunked replay archive:

```bash
.venv/bin/python -m microthermo precalculate \
  --preset carnot_discs --particles 200 --seed 123 --duration 0.1 \
  --fps 60 --chunk-frames 256 --output runs/disc_replay
```

The output directory must be empty. `ReplayReader` can load an exact saved
frame or seek to the latest frame at or before a physical time without
running the solver. To play it in the desktop view, run
`.venv/bin/python -m microthermo replay runs/disc_replay`. The
[format specification](docs/REPLAY_FORMAT.md) describes the saved columns,
interpolation limits, and short size and display probes.

The solver keeps swept-box candidates as sorted NumPy arrays, avoiding a
temporary Python tuple list. The public `swept_pairs_sweep` function still
returns its original list format.
With Numba and `--pair-kernel auto`, an exact swept-circumcircle check filters
pairs before polygon or disc CCD. `--pair-kernel scalar` retains the earlier
conservative reach filter for comparison.
For polygon pairs, the batch evaluates the separating-axis gap during most
advancement steps and computes full contact features when a pair reaches
contact. The benchmark warms the complete collision search before timing.

For a short disc throughput check with hardware and shaft speed recorded:

```bash
.venv/bin/python -m microthermo.measurements.performance_benchmark \
  --preset carnot_discs --particles 200 --duration 0.1 \
  --shaft-speed 0.15 --pair-search sweep --numeric-backend numba \
  --cam-kernel auto
```

To measure both gases at 16, 32, 64, 128, 200, and 500 particles, run the
isolated-process matrix. Each case saves its own warm solver time, event rate,
startup time, peak process memory, seed, radius, occupied area fraction,
backend, CPU, and numerical health. The report is updated after every case,
so completed measurements survive an interrupted run. The standard physical
durations shrink with particle count; compare throughput only alongside the
reported duration and event count.

```bash
.venv/bin/python -m microthermo.measurements.performance_matrix \
  --numeric-backend numba --pair-search sweep \
  --output runs/performance_matrix.json
```

Use `--presets carnot_triangles --counts 64,128,200 --seeds 123,124
--repeats 3` for a smaller repeated study. A nonstandard count requires
`--duration`. Each row runs in a fresh Python process, so process startup and
JIT compilation are visible separately from warm solver timing. These short
probes characterize throughput; they do not establish long-run collision
health or interactive playback.

For the longer three-seed 200/500 particle comparison used in the gap plan:

```bash
.venv/bin/python -m microthermo.measurements.performance_matrix \
  --counts 200,500 --seeds 123,124,125 --duration 0.1 \
  --numeric-backend numba --pair-search sweep \
  --output runs/performance_matrix_multiseed.json
```

For two complete cycles at the 96-triangle GUI default, with numerical
health recorded after every branch and a sampled/checkpoint comparison for
seed 123:

```bash
.venv/bin/python -m microthermo.measurements.triangle_cycle_health \
  --particles 96 --seeds 123 124 125 --cycles 2 \
  --shaft-speed 0.15 --verify-seed 123 --sample-cadence 0.83 \
  --output runs/triangle_default_96_two_cycles.json
```

This command requires Numba and saves each seed to the JSON report as it
finishes. Its sampling check compares exact events, body state, ledger, and
cycle markers; the longer run is separate from the short throughput matrix.

For one complete hot branch at 200 triangles across three seeds, use
`--branches 1` (each branch lasts `pi / (2 * shaft_speed)`):

```bash
.venv/bin/python -m microthermo.measurements.triangle_cycle_health \
  --particles 200 --seeds 123 124 125 --branches 1 \
  --shaft-speed 0.15 --verify-seed 123 --sample-cadence 0.83 \
  --output runs/triangle_200_hot_branch.json
```

For the larger 500-triangle full-branch health and throughput study, omit the
extra sampled/replay run with `--no-verify`:

```bash
.venv/bin/python -m microthermo.measurements.triangle_cycle_health \
  --particles 500 --seeds 123 124 125 --branches 1 \
  --shaft-speed 0.15 --no-verify \
  --output runs/triangle_500_hot_branch.json
```

For complete default-speed 200-triangle cycles, including a sampled replay
for seed 123 and Python/Numba next-collision comparisons from identical
checkpoints at the midpoint of each branch:

```bash
.venv/bin/python -m microthermo.measurements.triangle_cycle_health \
  --particles 200 --seeds 123 124 125 --cycles 1 \
  --shaft-speed 0.15 --verify-seed 123 --sample-cadence 0.83 \
  --same-state-queries --output runs/triangle_200_full_cycle.json
```

For one complete 500-triangle cycle with the four branch-midpoint backend
queries, use:

```bash
.venv/bin/python -m microthermo.measurements.triangle_cycle_health \
  --particles 500 --seeds 123 --cycles 1 --shaft-speed 0.15 \
  --no-verify --same-state-queries \
  --output runs/triangle_500_full_cycle_seed123.json
```

For the extended 200-triangle, two-cycle check across three seeds, including
a different sampling cadence and checkpoint replay for seed 124:

```bash
.venv/bin/python -m microthermo.measurements.triangle_cycle_health \
  --particles 200 --seeds 123 124 125 --cycles 2 \
  --shaft-speed 0.15 --verify-seed 124 --sample-cadence 1.37 \
  --same-state-queries \
  --output runs/triangle_200_two_cycles_3seed.json
```

To extend the 500-triangle full-cycle numerical and local backend check to a
second seed without the extra sampled trajectory:

```bash
.venv/bin/python -m microthermo.measurements.triangle_cycle_health \
  --particles 500 --seeds 124 --cycles 1 --shaft-speed 0.15 \
  --no-verify --same-state-queries \
  --output runs/triangle_500_full_cycle_seed124.json
```

For exact sampled/checkpoint replay validation over a complete 500-triangle
cycle, including a direct trajectory for comparison:

```bash
.venv/bin/python -m microthermo.measurements.triangle_cycle_health \
  --particles 500 --seeds 123 --cycles 1 --shaft-speed 0.15 \
  --verify-seed 123 --sample-cadence 0.83 \
  --output runs/triangle_500_full_cycle_sampled_seed123.json
```

This repeats the full solver run for the sampled trajectory and the latter
half after restoring a checkpoint, so allow substantially more time than the
direct-only 500-triangle command.

To check the first branch after one complete 500-triangle cycle, including
Python/Numba next-collision comparisons at all five branch midpoints:

```bash
.venv/bin/python -m microthermo.measurements.triangle_cycle_health \
  --particles 500 --seeds 123 --branches 5 --shaft-speed 0.15 \
  --no-verify --same-state-queries \
  --output runs/triangle_500_five_branches_seed123.json
```

On an i7-11850H this direct check took about 23 minutes; it does not repeat
the trajectory with sampling or checkpoint replay.

For a longer 500-triangle run, save the direct simulation at each branch end:

```bash
.venv/bin/python -m microthermo.measurements.triangle_cycle_health \
  --particles 500 --seeds 123 --cycles 2 --shaft-speed 0.15 \
  --no-verify --same-state-queries \
  --branch-checkpoint runs/triangle_500_two_cycles_seed123.branch.pkl.gz \
  --output runs/triangle_500_two_cycles_seed123.json
```

Run the same command again after an interruption to resume from the last
completed branch. The checkpoint uses Python pickle; load only files this
command created locally. Use a distinct checkpoint path for each seed and
configuration. The final JSON contains the completed measurements.

To repeat the two-cycle direct check for seed 124, change `--seeds 123` to
`--seeds 124` and use seed-124 paths for both `--branch-checkpoint` and
`--output`. The seed-124 run took about 36 minutes on an i7-11850H.

To compare a completed seed-124 direct checkpoint with a two-cycle sampled
trajectory at a different cadence and a midpoint replay:

```bash
.venv/bin/python -m microthermo.measurements.triangle_cycle_health \
  --particles 500 --seeds 124 --cycles 2 --shaft-speed 0.15 \
  --verify-seed 124 --sample-cadence 1.37 --same-state-queries \
  --branch-checkpoint runs/triangle_500_two_cycles_seed124.branch.pkl.gz \
  --output runs/triangle_500_two_cycles_sampled_seed124_cadence1p37.json
```

The checkpoint must have been made by the direct command above with the
same seed and configuration. This skips the direct solve but still computes
the full sampled trajectory and the post-midpoint replay. On an i7-11850H,
these took about 38 and 22 minutes, respectively, for two cycles.

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
