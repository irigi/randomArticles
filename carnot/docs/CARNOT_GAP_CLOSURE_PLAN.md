# Carnot gap closure: remaining work

Updated: 2026-10-05

The gap is not closed. The application has a working experimental apparatus,
plots, diagnostics, event accounting, and a reproducible speed study. The full
change history and completed tests are in `docs/PROGRESS.md`; this file tracks
open gates and their measurement evidence. The 24-run baseline and 16-run
optional cold-jacket studies support a hot-heat speed effect for 16 discs and
improved cold-sector coupling for both gas types, but all drift screens are
inconclusive and no stationary efficiency trend is established.
The measured numerical and desktop gates now pass, while the stationary
Carnot thermodynamic claim and release acceptance remain open.

## Remaining external run-and-review iterations (updated after dilution pilot)

| Work | Next iterations | Condition |
| --- | ---: | --- |
| Physical protocol diagnosis and redesign | No external batch queued | Dilution failed its fixed rule; use saved evidence before another run. |
| New thermodynamic pilot and possible extension/comparison | At least 1, conditionally 2+ | A justified protocol must be specified first; no finite upper bound is known. |
| Clean-OS release check | 1 | Bash launcher prepared for a separate clean Ubuntu 24.04 checkout; full brief map also remains. |

The dilution batch is complete and is no longer counted. There are **zero
science runs currently queued**. From this point, at least **two external
iterations** remain on a successful path: one new thermodynamic pilot and one
clean-OS release check. A favorable pilot would likely add a longer
stationarity/comparison batch, making **three or more**. The physics redesign
has no reliable upper bound; this estimate will change only when a new
protocol and its decision rule are justified. In-session analysis, code work,
and requirement mapping are separate from these user-run batches.

## Next checkpoints

These are the actionable steps within the larger workstreams below. Completed
measurements remain visible here so progress does not look like a repeating
five-item queue.

- [x] Three-seed, two-cycle numerical health at the 96-triangle default.
- [x] Three-seed, complete hot-branch numerical health at 200 triangles,
  including one sampled run and checkpoint replay.
- [x] Run a complete hot branch at 500 triangles across three seeds, recording
  throughput, collision health, radius, and occupied fraction.
- [x] Run complete default-speed cycles at 200 triangles across seeds, with
  sampled/checkpoint replay and same-state backend queries.
- [x] Extend 500-triangle validation beyond the hot branch and compare
  same-state backend queries; document the resource cost of full cycles.
- [x] Profile and improve the measured high-count bottlenecks; report whether
  200 triangles reach interactive playback and whether 500 remains below it.
- [x] Extend 200-triangle numerical health to two complete cycles across three
  seeds, with branch-midpoint same-state backend queries and a second sampled
  seed/cadence plus checkpoint replay.
- [x] Extend 500-triangle full-cycle numerical health and branch-midpoint
  same-state backend queries to a second seed.
- [x] Validate a complete 500-triangle cycle with regular sampling and
  checkpoint restore against an exact direct event history and ledger.
- [x] Extend seed-123 500-triangle collision health through the first branch
  after a complete cycle, with five branch-midpoint backend queries.
- [x] Complete two 500-triangle cycles for seed 123 with branch-end recovery
  and eight branch-midpoint backend queries.
- [x] Extend two-cycle 500-triangle collision health and branch-midpoint
  backend queries to a second seed.
- [x] Check two-cycle 500-triangle sampling and checkpoint replay at another
  cadence against the exact direct event history and ledger.
- [x] Complete the 200/500-triangle baseline measurements. Further unchanged
  long runs are unnecessary: the tested cases pass collision-health and
  sampling checks, while live throughput misses the section 1 target.
- [x] Defer further live-solver optimization at the user's request. Record
  the missed live 200-triangle target and pursue precalculated replay for
  smooth high-count display; recheck parity and throughput only after a
  future solver change.
- [x] Collect and interpret multi-seed, 18-cycle thermodynamic evidence for
  discs and triangles at multiple shaft speeds.
- [x] Diagnose overheated hot-branch entry and test an optional cold-sector
  jacket in matched 18-cycle disc and triangle ensembles, keeping the original
  apparatus and screening gates as controls.
- [x] Extend the slow 32-triangle cold-jacket cohort to 36 cycles across three
  seeds without changing the heat or drift gates; distinguish positive hot
  input from stationary behavior.
- [x] Test increased hot-sector contact area against the matched slow
  16-disc cold-jacket control. Contacts increased about fourfold, but paired
  hot-heat changes had mixed signs and all drift screens stayed inconclusive.
- [x] Test longer compression and shorter expansion against the matched slow
  16-disc cold-jacket control, preserving hot/cold exposure lengths. Paired
  hot-heat changes had mixed signs, heat-gate readiness fell to 1/3, and all
  drift screens stayed inconclusive.
- [x] Test longer cold exposure with shorter adiabatic expansion against the
  matched slow 16-disc cold-jacket control. Hot heat fell in all three seeds,
  heat-gate readiness fell to 1/3, and all drift screens stayed inconclusive.
- [x] Compare the existing 16- and 32-disc slow cold-jacket artifacts at
  matched seeds and screens. Heat-gate readiness improved from 2/3 to 3/3;
  per-particle hot-heat changes had mixed signs and drift remained inconclusive.
- [x] Extend two-cycle 96-triangle numerical health to seeds 126–127 at two
  sampling cadences, including checkpoint replay and same-state backend
  queries; both passed.
- [x] Test hot-plus-cold jackets for 32 discs at two shaft speeds and matched
  seeds. All six runs passed the positive-heat block gate, but all drift
  screens remained inconclusive and the earlier speed ordering was not robust.
- [x] Precalculate ten-physical-second 200/500-disc replay archives and
  measure storage and offline write time. Real-display acceptance remains open.
- [x] Probe forward and reversed free-shaft trajectories at two seeds; none
  completed a cycle, so the cycle-level comparison remains open.
- [x] Extend matched both-jacket 32-disc cohorts to 36 cycles at two speeds.
  All six runs retained positive-hot-heat readiness, while all twelve drift
  screens remained inconclusive; 144 unit tests and scientific checks passed.
- [x] Diagnose the saved 36-cycle branch records and predeclare a matched
  compression-timing pilot with fixed heat and drift screens. See
  `docs/CARNOT_NEXT_PROTOCOL_2026-10-04.md`.
- [x] Run the predeclared 18-cycle, two-speed, three-seed compression-timing
  pilot. It missed the hot-entry-temperature decision threshold at both
  speeds, so stop this variant without a longer run.
- [x] Run the predeclared 64-disc, two-speed, three-seed both-jacket pilot.
  All six runs passed numerical and positive-heat checks, but zero of twelve
  drift screens were bounded; stop this variant without a longer extension.
- [x] Run the predeclared 32-disc dilution pilot at two radius scales, two
  speeds, and three seeds. The quarter-scale hot-entry threshold failed at
  fast speed; no drift screen was bounded. Stop this variant.
- [x] Benchmark sustained real-display playback and seek/load behavior for the
  saved 200/500-disc replay archives. Both averaged about 61 paint/s, but
  three pauses above 25 ms per archive leave smooth playback open.
- [x] Remove long paint gaps and reduce seeks for the ten-second 200/500-disc
  archives with the warmed linear reader: both rendered at 62.5 paints/s
  with no gap above 25 ms in that check.
- [x] Reconstruct disc positions and velocities between frames from recorded
  impulses. At 200/500 discs, the one-second midpoint errors fell to roundoff
  level and all 146 tests passed.
- [x] Diagnose paint gaps, stress 19-chunk seeks, and test 1×/4× display at
  both disc counts. Random uncached seeks stayed below 10 ms; 1× garbage
  collections coincided with the long paint gaps. Triangle midpoint error
  was material in archives without contact points.
- [x] Verify playback garbage-collection control and event-aware triangle
  reconstruction at 32/96 triangles, seed 123. All six desktop cases had no
  paint gap above 25 ms; 120 triangle midpoints per count matched direct
  simulation to roundoff. The 147-test suite and scientific checks passed.
- [x] Extend triangle replay to seed 124 for five physical seconds at 32/96
  particles; compare 300 midpoints, saved statistics, and ledgers with the
  direct source, then check desktop playback and UI seeks.
- [x] Verify streaming event export's high-count memory and exact archive
  parity, then measure genuine 19-chunk desktop playback. All checks passed
  at 200/500 discs; the full 148-test suite and scientific checks passed.
- [x] Bound default controlled-cam apparatus interpolation error analytically
  and make selector-shoe jumps occur at logged branch transitions. Focused
  replay tests pass; see `docs/REPLAY_FORMAT.md`.
- [x] Verify 200/500-triangle archive fidelity through a branch transition:
  source motion, saved statistics, branch/selector state, and apparatus
  geometry pass; the 149-test suite and scientific checks pass.
- [x] Diagnose first-paint and play-to-end timing for high-count triangle 4×
  desktop playback. Four repeated 4× runs reached the end at the expected
  wall time with 166 ticks and paints and no gap above 17.5 ms.
- [ ] Establish stationary positive hot input across seeds and shaft speeds;
  then compare efficiency, pressure–area work, and free-shaft work against the
  reference under the gate in section 2.
- [x] Finish numerical acceptance under the measured gate in section 3:
  long dilute runs, multi-feature clusters, moving walls, rotating-pair CCD,
  tolerance convergence, and sampling/checkpoint parity.
- [x] Add desktop seed reset, fixed-duration and branch steps, target playback
  rate, recorded physical shaft-speed intervention, preset switching, and
  zoom/pan/particle selection. Focused UI and intervention tests pass.
- [x] Verify the first control pass on the real desktop against an identical
  headless command schedule.
- [x] Add configuration, checkpoint, run and plot file actions, a refining
  status, and a two-row toolbar; the focused offscreen round trip passes.
- [x] Complete the measured desktop controls and presentation gate in section 4:
  headless parity, file actions, real-display response, recording, inspector,
  and dropped-render-frame physics parity.
- [ ] Complete release acceptance under the gate in section 5.
- [ ] Run the clean Ubuntu installation check from
  `docs/CARNOT_CLEAN_UBUNTU_PROTOCOL_2026-10-05.md` and review both isolated
  install paths before deleting its temporary batch.
- [x] Add precalculated smooth replay under the measured gate in section 6
  for controlled-shaft Carnot discs and triangles at 200/500 particles.

## 1. Make larger particle runs practical

Selected route: **A+B**, conservative spatial search plus optional Numba
kernels. The implemented stages keep a fixed,
dilute occupied area fraction as count changes (triangles 3.12%, discs 7.54%
of the minimum cylinder area); a swept circumcircle grid replaces the dense
pair search from 32 particles; the numerical reach filter and full polygon
pair conservative advancement have Numba paths; and cam branch boundaries
are cached. Swept circumcircle bounds prune wall queries, and Numba kernels
now handle fixed-velocity polygon wall CCD and batched wall penetration
bounds. A selectable compiled swept-box search reduces candidate construction
cost at 200/500 particles. Triangle pair CCD can now be batched across
candidate pairs. A compiled wall reach mask and batched disc CCD now reduce
larger disc runs. The grid, scalar pair query, unpruned wall, and
all-pairs/Python paths remain selectable for
reference comparisons. This is a measured partial speedup,
not completion of the 200/500-particle performance gate.

- Extend performance and collision-health measurements to full branches and
  cycles once backend parity is established. The 16–500 particle matrix now
  covers three seeds, with a 0.1-physical-second window at 200/500 particles;
  full hot branches are now measured at 200 and 500 triangles across three
  seeds. Complete two-cycle studies now cover three 200-triangle seeds and two
  500-triangle seeds, with sampled/checkpoint comparisons at both counts.
- Profile the remaining scheduler and geometry costs at 200/500 particles.
  The fixed-velocity wall loop and wall penetration bounds are compiled;
  the cam piston still uses Python motion queries. A short 200-triangle
  profile now spreads time across pair-call orchestration, swept-grid
  construction, wall dispatch, and penetration. The swept-box search and
  batch triangle pair queries are measured improvements. Larger disc profiles
  show repeated cam piston motion queries and pair calls as current hot spots;
  the compiled wall reach mask and disc pair batch address part of this.
  Controlled-speed cam piston queries now use a specialized compiled batch.
  Fixed-wall disc queries and disc-only penetration now run in compiled
  batches. Swept-pair candidates now stay in arrays through the solver.
  Batched disc pairs now decode only contacts near the next event. Scheduler
  checkpoints share immutable cam geometry and copy ledger counters by value;
  no-portal energy checks skip region accounting. Triangle overlap diagnostics
  now batch compiled SAT gaps. Controlled-speed triangle contacts with the
  cam piston now use a compiled batch. An exact swept-circumcircle filter
  reduces polygon CCD calls before the batch. Polygon CCD now uses a SAT-only
  gap during most advancement steps and defers contact decoding until the next
  event is known. Packed triangle geometry now avoids repeated stacking and
  reduces checkpoint copy overhead. Triangle overlap diagnostics now use the
  SAT gap without decoding contact witnesses, and polygon pair CCD reuses
  vertex scratch arrays across advancement steps. Triangle-only penetration
  now scans sorted pair candidates and infinite walls in one compiled kernel;
  benchmark startup includes its JIT warmup. Constant-velocity triangle wall
  queries now run in a compiled body/wall batch with the scalar conservative
  bound preserved. Full-hot-branch backend parity now passes at 12 triangles
  for two seeds and two sampling cadences, including checkpoint replay.
  Same-state queries at three later points also match; quantify numerical
  sensitivity and extend long-run parity using collision-health and
  observable tolerances, alongside local event parity checks. Same-state
  contact queries now match in all four branches for two 32-triangle seeds;
  their full-cycle compiled runs are sampling-independent. A three-seed
  32-triangle full-cycle check now passes at the default shaft speed. Extend
  collision-health and parity acceptance beyond the current two-cycle
  200/500-triangle results only when a changed solver or a new failure case
  warrants it.
  Keep `fastmath` disabled and the Python path.
- Improve particle placement at 200/500 particles if setup time warrants it.
  Confirm the radius rule, density, and finite-size effects in every report.
- Keep the swept-grid bounds conservative for rotating bodies; add targeted
  high-speed, grazing, boundary, and moving-piston parity cases. Extend
  bounded-wall parity across longer moving-piston trajectories and seeds.
- At every stage compare same-state backend collision queries and event
  chronology while trajectories remain close, then compare heat/work
  observables, penetration, checkpoint/restart, and failure diagnostics
  against the Python reference across fixed seeds, triangle edge cases, and
  multiple sample cadences. Measure speed only after the physics gate passes.

Gate: the chosen backend provides a measured improvement on representative
Carnot runs, with unchanged scientific results within declared tolerances. Aim
for 200 triangles at an interactive achieved playback rate, then validate 500
at a documented density and hardware configuration. Publish the benchmark and
any cases that miss the target. Compilation time and fallback behavior count
as part of the result.

The user chose to defer the unmet live-rate target after the completed
200/500-triangle baseline studies. It is recorded as an unmet target, not a
passed gate; no further unchanged high-count performance runs are planned.
Smooth high-count display is now assigned to the precalculated replay path in
section 6.

The 2026-09-30 isolated-process matrix is saved in
`runs/performance_matrix_2026-09-30.json`. It covers both presets at 16, 32,
64, 128, 200, and 500 particles, seed 123, shaft speed 0.15, swept-pair search,
and Numba on an i7-11850H. Each row records warm solver time, throughput,
events, setup and JIT startup, process time, peak RSS, radius, occupied area
fraction, penetration, failures, and energy residual. Warm achieved physical
seconds per wall second for discs were 33.3, 15.6, 4.87, 1.94, 0.836, and
0.112, respectively; for triangles 28.1, 8.77, 2.42, 0.519, 0.173, and
0.0375. The 200-triangle probe advanced 0.02 physical seconds in 0.116 wall
seconds (43 events), and 500 triangles advanced 0.005 in 0.133 wall seconds
(14 events). All rows reported zero CCD failures and negligible penetration;
the windows are too short to establish long-run health. The matrix runner
records each case in a fresh process and saves after each case.

The three-seed extension is in
`runs/performance_matrix_multiseed_small_2026-09-30.json` (16–128 particles,
standard durations) and `runs/performance_matrix_multiseed_2026-09-30.json`
(200/500 particles, 0.1 physical seconds per run). At 200 particles, median
warm physical/wall throughput was 0.822 for discs and 0.191 for triangles;
at 500 it was 0.123 and 0.0160. Triangle ranges across seeds were
0.153–0.195 at 200 and 0.0157–0.0176 at 500. All 36 rows reported zero CCD
failures, maximum reported penetration below `1.9e-16`, and absolute energy
residual below `4.4e-13`. These short, forward-only trajectories support
throughput comparison but do not establish long-run collision health or
backend parity. The 500-triangle 0.1-second runs were about 60 times slower
than real-time playback on the measured machine.

A warmed `cProfile` run over 0.05 physical seconds assigned 1.52 of 3.55
profiled wall seconds to polygon pair CCD and 0.93 to penetration at 500
triangles; 1,501 repeated polygon stacks cost 0.28 seconds cumulatively.
Triangle vertices are now packed once into an `(n, 3, 2)` array, and
checkpoint copies clone that array in one operation while keeping each copy's
geometry independent. A comparable profile took 3.22 seconds overall, with
polygon pair CCD at 1.42 and penetration at 0.87 seconds. On the same
three-seed 0.1-second probes, median achieved physical/wall throughput rose
from 0.191 to 0.212 at 200 triangles and from 0.0160 to 0.0180 at 500.
The post-change measurements are in
`runs/performance_matrix_packed_triangles_2026-09-30.json`. Event counts,
energy residuals, penetration, and CCD failures matched the saved baseline
for all six cases. Polygon CCD and Python wall/penetration dispatch remain
large costs; interactive playback is still far away.

The next geometry pass removed repeat vertex-array allocations inside each
polygon pair CCD query and skipped witness decoding in triangle overlap
diagnostics. On the same three-seed, 0.1-physical-second triangle probes,
median achieved physical/wall throughput increased from 0.212 to 0.242 at
200 particles and from 0.0180 to 0.0226 at 500. Event counts, maximum
penetration, CCD failure counts, and energy residuals matched all six saved
baseline cases. Results are in
`runs/performance_matrix_scratch_vertices_2026-09-30.json`; an intermediate
SAT-only diagnostic measurement is in
`runs/performance_matrix_sat_gap_2026-09-30.json`. The 500-particle runs
remain about 44 times slower than real time, and full-branch parity is open.

A subsequent warmed 500-triangle profile over 0.03 physical seconds assigned
0.44 of 1.59 profiled wall seconds to penetration, 0.62 to polygon pair CCD,
and 0.15 to wall dispatch. Fusing the body-major triangle pair and infinite-wall
penetration scan into one Numba kernel reduced penetration to 0.11 and total
profiled time to 1.24 seconds on the same seed and interval. The six-case
three-seed benchmark shows median physical/wall throughput rising from 0.242
to 0.294 at 200 triangles and from 0.0226 to 0.0270 at 500, with identical
event counts, energy residuals, maximum penetration, and CCD failure counts.
The artifact is `runs/performance_matrix_triangle_penetration_2026-09-30.json`.
The 500-triangle window is still about 37 times slower than real time;
full-branch parity and the interactive target remain open.

The next wall-dispatch pass batched triangle CCD against constant-velocity
infinite walls, retaining the scalar wall reach bound, contact features,
failure reasons, and the selectable Python path. On the same three-seed
0.1-second probes, median physical/wall throughput rose from 0.294 to 0.349
at 200 triangles and from 0.0270 to 0.0300 at 500. Event counts, energy
residuals, maximum penetration, and CCD failures matched all six preceding
runs. The artifact is
`runs/performance_matrix_triangle_wall_batch_2026-09-30.json`. A comparable
500-triangle profile over 0.03 physical seconds fell from 1.24 to 1.14
profiled seconds; the 7,713 per-body `_wall_toi` calls disappeared from the
hot list. Polygon pair CCD now dominates at 0.61 seconds. The 500-triangle
short window remains about 33 times slower than real time; full-branch parity
and the interactive target remain open.

Full-hot-branch triangle parity now has a regression test at 12 particles,
seeds 123 and 124, prescribed shaft speed 1.5, and two sampling cadences.
It compares ordered event kinds/participants/times, contact-related event
metadata and impulses, final body state, every ledger counter, penetration,
energy residual, failure diagnostics, and checkpoint/restart. Both seeds
match within declared tolerances and have zero CCD failures and zero reported
penetration. The accelerated shaft makes the Python reference practical; it
does not establish default-speed or 200/500-particle parity. A separate
full-cycle pilot at the same settings diverged in event chronology after the
first branch: seed 123 first differs around event 65 (about 2 physical
seconds), and seed 124 around event 78. Both finish without CCD failures,
  and cannot be treated as full-cycle event parity. The first material timing
  drift for seed 123 appears at event 42 (about 1.253 physical seconds), where
  event time differs by 4.6e-8 seconds. Independent states differ by about
  4.0e-6 in position/velocity at 1.24 seconds, 8.6e-3 at 1.75 seconds, and
  0.657 at 1.96 seconds. Re-querying from the exact same compiled checkpoint
  at those three times gives the same next contact, time, features, point, and
  normal in both backends. At 1.96 seconds, the independent Python trajectory
  predicts a different pair contact while Python restarted from the compiled
  state predicts the same wall contact as Numba. This supports amplified
  trajectory sensitivity as the cause of the observed later discrepancy;
  it does not prove every collision path equivalent. Extend same-state local
  parity cases and assess longer-run physical observables and numerical health
  across seeds; exact independent event chronology should be required only
  while trajectories remain within a declared state tolerance.

The next full-cycle regression uses 32 triangles at prescribed shaft speed
1.5 for seeds 123 and 124. At the midpoint of each of the four branches,
Python and Numba receive identical checkpoints and agree on the next contact
time, participants, features, point, and normal. These eight queries include
pair and fixed/moving-wall contacts. The compiled full cycles contain 664 and
690 events, respectively, with zero CCD failures, zero reported penetration,
and absolute energy residuals below `3.5e-13`. Direct and sampled runs at
cadences 0.083 and 0.137 have exactly equal events, final body arrays, ledger,
and cycle markers; a mid-cycle checkpoint reproduces the final result. This
is a stronger local/sampling gate, not an independent-trajectory event-parity
claim or a default-speed/high-count acceptance run.

At the default shaft speed 0.15, three 32-triangle, one-cycle runs (seeds
123–125) recorded 5,121, 5,584, and 4,925 events. Each had zero CCD failures,
zero reported penetration, and absolute energy residual below `1.5e-13`.
Python and Numba agreed on same-state next-collision queries at each branch
midpoint (12 queries total, including one shared no-collision result). Three
different sampling cadences and mid-cycle checkpoint replay reproduced the
compiled event histories, body states, ledgers, and cycle markers exactly.
This covers 32 particles for one cycle; larger and longer acceptance runs
are measured separately below.

The first longer 96-triangle default-speed study is saved in
`runs/triangle_default_96_two_cycles_2026-09-30.json`. Seeds 123–125 each
completed two cycles, with 44,607, 45,160, and 44,215 events; zero CCD
failures and zero reported penetration; and absolute energy residuals below
`2.4e-12`. Radius was 0.01768 and occupied area fraction 3.12% in every
run. Warm solver time was 36.4–37.2 wall seconds per 83.8 physical seconds
on an i7-11850H, or 2.25–2.30 physical seconds per wall second. Seed 123's
0.83-second sampling cadence and mid-run checkpoint replay matched its direct
events, body arrays, ledger, and cycle markers exactly. The new
`triangle_cycle_health` command saves each seed as it completes. This is
useful two-cycle collision-health evidence at the GUI default, not a
many-cycle thermodynamic or high-count performance gate. Sampling/replay was
checked for one seed; same-state Python parity at 96 and larger-count
long-run health remain open.

The next study, `runs/triangle_200_hot_branch_2026-09-30.json`, completed
the entire default-speed hot branch (10.47 physical seconds) for seeds
123–125 at 200 triangles. It recorded 22,484, 22,738, and 22,879 events,
zero CCD failures and zero reported penetration, and absolute energy
residuals below `1.5e-12`. Radius was 0.01225 and occupied area fraction
3.12%. Warm solver time was 35.0–36.6 wall seconds, or 0.286–0.299 physical
seconds per wall second on the i7-11850H. Seed 123's sampled run and
checkpoint replay matched its direct events, body state, ledger, and cycle
markers exactly. This closes the 200-triangle full-branch measurement
checkpoint, but not the interactive speed or full-cycle acceptance gates.

The three-seed 500-triangle full hot-branch measurement is saved in
`runs/triangle_500_hot_branch_2026-09-30.json`. At the default shaft speed,
seeds 123–125 completed 10.47 physical seconds with 77,834, 79,270, and
77,791 events. All reported zero CCD failures and zero penetration; absolute
energy residual stayed below `5.4e-12`. Radius was 0.007746 and occupied area
fraction was 3.12%. Warm solver time was 334.8–359.9 wall seconds, or
0.0291–0.0313 physical seconds per wall second on the i7-11850H. This
closes the three-seed 500-triangle full-branch measurement checkpoint. It
does not establish full-cycle health or parity across those three seeds.
Sampling/replay was not repeated at 500. The single-seed full-cycle result
below replaces the earlier first-branch runtime extrapolation.

The three-seed, one-cycle 200-triangle study is in
`runs/triangle_200_full_cycle_2026-09-30.json`. At the default shaft speed,
seeds 123–125 completed 41.89 physical seconds and 56,513, 59,324, and
58,209 events. Each run reported zero CCD failures and zero penetration;
absolute energy residual was below `2.4e-12` at every branch end. The radius
was 0.01225 and occupied area fraction 3.12%. Warm direct-solver time was
80.5–85.9 wall seconds, yielding 0.488–0.520 physical seconds per wall
second. All 12 Python/Numba next-collision queries from identical branch
midpoint checkpoints matched status, timing, contact count/features, point,
and normal within the declared tolerances. Seed 123's 0.83-second sampling
cadence and checkpoint replay exactly matched direct events, body state,
ledger, and cycle markers. Query time was excluded from the solver timing.
This closes the 200-triangle full-cycle health and local-parity checkpoint,
not the interactive speed or many-cycle thermodynamic gates.

The two-cycle extension is in
`runs/triangle_200_two_cycles_3seed_2026-09-30.json`. Seeds 123–125 completed
83.78 physical seconds with 117,550, 120,041, and 117,783 events. All three
reported zero CCD failures and zero penetration, with maximum absolute
branch-end first-law residual below `3.9e-12`. All 24 Python/Numba
next-collision queries from identical branch-midpoint states matched; the
largest time difference was below `2.7e-16` seconds. Seed 124's 1.37-second
sampling cadence and mid-run checkpoint replay exactly reproduced direct
events, body arrays, ledger, and cycle markers. Together with the earlier
seed-123 0.83-second one-cycle check, sampling now covers two seeds and
cadences at 200 triangles, though only one seed/cadence was repeated over two
cycles. Direct solver throughput was 0.424–0.440 physical seconds per wall
second on the i7-11850H. These are numerical-health and local backend-query
results; independent Python/Numba trajectories need not retain identical
event histories over chaotic multi-cycle evolution. Live throughput remains
below real time.

The seed-123 500-triangle full cycle is in
`runs/triangle_500_full_cycle_seed123_2026-09-30.json`. Its four branches
completed 41.89 physical seconds and 198,943 events with zero CCD failures
and zero reported penetration. The maximum absolute branch-end energy
residual was `1.36e-11`. The radius was 0.007746 and occupied area fraction
3.12%. Direct solver time was 823.6 wall seconds (13.7 minutes), excluding
the four local backend queries, for 0.0509 physical seconds per wall second
on the i7-11850H. All four Python/Numba next-collision queries from
identical branch-midpoint checkpoints matched status, timing, contacts,
features, point, and normal within the declared tolerances. The first branch
reproduced the earlier seed-123 500-triangle hot-branch event count and
energy residual.

The second full-cycle 500-triangle seed is in
`runs/triangle_500_full_cycle_seed124_2026-09-30.json`. Seed 124 completed
41.89 physical seconds with 205,758 events, zero CCD failures, zero reported
penetration, and a maximum absolute branch-end first-law residual of
`1.30e-11`. All four branch-midpoint Python/Numba next-collision queries from
identical states matched. Its direct solver took 924.7 wall seconds, or
0.0453 physical seconds per wall second, compared with 0.0509 for seed 123.
Both runs used radius 0.007746 and 3.12% occupied area fraction. The two
full-cycle seeds strengthen numerical-health and local-query evidence, but
500-triangle sampling, multiple cycles, and interactive live throughput
remained open at that checkpoint. The sampled result below closes the
single-cycle sampling check.

The seed-123 500-triangle sampling study is in
`runs/triangle_500_full_cycle_sampled_seed123_2026-09-30.json`. At a 0.83
physical-second sample cadence, its complete-cycle sampled run and a replay
from a checkpoint just after the midpoint reproduced the direct event list,
body arrays, ledger, and cycle markers exactly. The new direct trajectory
also reproduced the earlier seed-123 artifact's 198,943 events, event-kind
counts, branch records, ledgers, and residuals exactly. It had zero CCD
failures, zero reported penetration, and maximum absolute branch-end
first-law residual `1.35e-11`. Direct, sampled, and restored solver times
were 859.9, 799.6, and 346.6 wall seconds, respectively, on the i7-11850H.
This is a numerical checkpoint/restart test, not smooth visual replay.
Further 500-triangle multi-cycle and multi-cadence acceptance remains open;
the direct rate of 0.0487 physical seconds per wall second remains far below
real-time live playback.

The first 500-triangle post-cycle branch is in
`runs/triangle_500_five_branches_seed123_2026-10-01.json`. Seed 123 completed
five branches (52.36 physical seconds) with 288,093 events, zero CCD failures,
zero reported penetration, and a maximum absolute branch-end first-law
residual of `1.35e-11`. The first four branch records exactly reproduce the
earlier one-cycle seed-123 result (198,943 events). All five same-state
Python/Numba next-collision queries matched. The direct solver took 1,381.1
wall seconds on the i7-11850H, or 0.0379 physical seconds per wall second;
radius and occupied fraction remained 0.007746 and 3.12%. This extends
collision-health evidence into the next hot branch, but it is not a second
complete cycle or a sampled/checkpoint comparison beyond the first cycle.

The seed-123 500-triangle two-cycle result is in
`runs/triangle_500_two_cycles_seed123_2026-10-01.json`. Eight complete
branches covered 83.78 physical seconds and 416,504 events. The first four
and first five branch records exactly reproduced the earlier one-cycle and
five-branch artifacts, respectively. All eight branch-midpoint same-state
Python/Numba collision queries matched. There were zero CCD failures, zero
reported penetration, and a maximum absolute branch-end first-law residual
of `1.35e-11`. Direct solver time was 2,122.2 wall seconds on the i7-11850H,
or 0.0395 physical seconds per wall second. Radius remained 0.007746 and
minimum-area occupancy 3.12%. The runner now saves a compressed branch-end
state for interruption recovery; the final file was 17.6 MB. A focused
interruption/resume test reproduced the uninterrupted event, ledger, and
diagnostic results. This one-seed direct run does not establish multi-seed,
multi-cadence sampling independence or interactive throughput.

The second 500-triangle two-cycle seed is in
`runs/triangle_500_two_cycles_seed124_2026-10-01.json`. Seed 124 completed
eight branches (83.78 physical seconds) with 421,832 events. The first four
branch records exactly reproduced its earlier one-cycle artifact. All eight
branch-midpoint same-state Python/Numba collision queries matched. There were
zero CCD failures, zero reported penetration, and a maximum absolute
branch-end first-law residual of `1.30e-11`. Direct solver time was 2,169.5
wall seconds on the i7-11850H, or 0.0386 physical seconds per wall second.
The radius and minimum-area occupancy again were 0.007746 and 3.12%. Together
with seed 123, this gives two-seed, two-cycle direct collision-health evidence;
sampling independence beyond one cycle and interactive throughput remain open.

The two-cycle seed-124 sampling study is in
`runs/triangle_500_two_cycles_sampled_seed124_cadence1p37_2026-10-01.json`.
It reused the completed direct branch checkpoint, then advanced a separate
500-triangle trajectory at a 1.37-physical-second sample cadence and replayed
from a checkpoint just after the midpoint. Both paths exactly matched the
direct 421,832-event history, body arrays, ledger, and cycle markers. The
restored direct result also exactly matched the earlier seed-124 two-cycle
artifact, including all eight branch records and same-state backend queries.
The sampled solve took 2,268.4 wall seconds and the post-midpoint replay
1,336.5 seconds, versus 2,169.5 seconds for the direct solve on the
i7-11850H. This extends sampling/checkpoint evidence from one cycle at 0.83
seconds to two cycles at 1.37 seconds. It does not provide smooth visual
replay or interactive live performance.

Baseline high-count measurement is complete for the tested default-speed
configuration. The three 200-triangle two-cycle seeds achieved 0.424–0.440
physical seconds per wall second; the two 500-triangle two-cycle seeds
achieved 0.0386–0.0395. Repeating these unchanged long runs is unlikely to
alter the performance conclusion. The section 1 live target remains unmet;
future measurements should be tied to a specific solver change, new physical
configuration, or regression, with precalculated replay tracked separately.

A fresh 0.1-physical-second profile assigned 2.26 of 4.30 profiled wall
seconds to compiled polygon pair CCD at 500 triangles, versus 0.29 to swept
candidate construction and 0.45 to penetration checks. The pair batch now
reuses two triangle-vertex scratch arrays across candidates, preserving the
scalar CCD path. The isolated-process before/after matrices are
`runs/performance_matrix_pair_scratch_before_2026-09-30.json` and
`runs/performance_matrix_pair_scratch_after_2026-09-30.json`, each with three
seeds at 200 and 500 triangles for 0.1 physical seconds. Median warm
physical/wall throughput at 500 rose from 0.0285 to 0.0300 (+5.4%); at
200 it changed from 0.335 to 0.331 (-1.2%, within timing variation).
Event counts, event-kind counts, ledgers, CCD failures/refinements, maximum
penetration, and energy residuals matched exactly in all six cases. Median
JIT warmup remained about 0.13 seconds with cached compiled kernels; median
process startup was 0.17 seconds at 200 and 0.38 at 500. The Python backend
remains available. All 129 tests pass.

The 200-triangle direct solver is still below interactive real time: the
measured full hot branches achieved 0.286–0.299 physical seconds per wall
second and full cycles 0.488–0.520. The one 500-triangle full cycle achieved
0.0509. This profiling checkpoint is complete, but the section 1 interactive
target is missed. The precalculated replay work in section 6 is the planned
route to smooth display at these particle counts; any additional live-solver
speedups must pass the same numerical parity gates.

First-slice measurements (i7-11850H, Python 3.12.3, NumPy 2.4.6, Numba
0.67.0, seed 123, new radius rule, warm JIT, short physical intervals):
32 discs/0.3 s: all-pairs Python 0.388 s, grid Python 0.267 s, grid Numba
0.264 s; 64 discs/0.15 s: 1.584, 0.798, 0.784 s. Before compiling polygon
witnesses, 32 triangles/0.15 s took 2.453 s on grid Python and 2.452 s on
grid Numba; 64 triangles/0.08 s took 8.504 and 8.494 s. With the polygon
witness kernel compiled and warmed, 32 triangles/0.15 s took 2.493 s on grid
Python and 0.314 s on grid Numba; 64 triangles/0.08 s took 8.871 and 1.038 s.
Event counts matched within each row. These short samples are directional,
not long-run throughput evidence. The compiled witness kernel yields roughly
8× in these short triangle runs; conservative advancement and wall handling
still need profiling at larger counts.
Triangle gases at 200 and 500 particles initialize without overlap under the
new radius rule (about 0.05 s and 0.28 s setup on the same machine); their
long-run collision health remains unmeasured. Fusing polygon transforms into
the compiled witness kernel reduced a 200-triangle, 0.02-physical-second run
from 7.06 to 3.72 wall seconds (43 events). A 500-triangle,
0.005-physical-second run took 5.08 wall seconds (14 events). Neither is near
interactive playback. Both short runs had zero reported penetration and CCD
failures, with energy residuals below `1.4e-13`.

The next A+B slice compiled the complete polygon pair conservative-advancement
loop, including touching-contact re-query and impact-time bisection. On the
same short probes, 200 triangles took 2.38 wall seconds for 0.02 physical
seconds (43 events) and 500 triangles took 2.03 wall seconds for 0.005 physical
seconds (14 events), with zero reported penetration/CCD failures and residuals
below `1.4e-13`. Direct CCD edge cases and a sampled checkpoint/replay case
were added to the backend parity checks. These are still far from interactive
playback and do not establish long-run collision health.

The bounded-wall slice caches wall speed/position per interval and rejects a
wall only if the particle's swept circumcircle cannot reach its infinite
line. Penetration checks use circle lower bounds before polygon geometry;
the unpruned path remains selectable as `--wall-search all`. On the same
machine and seeds, 200 triangles/0.02 physical seconds took 2.18 wall
seconds unpruned and 1.00 bounded (43 events); 500 triangles/0.005 physical
seconds took 1.99 and 1.05 (14 events). Event counts, penetration, and energy
residuals matched. A short post-change profile at 200 triangles assigned
about 0.026 of 0.087 profiled seconds to the near-wall polygon loop, 0.020
to penetration, and 0.010 to grid construction. Longer parity and scientific
acceptance runs remain open.

The fixed-velocity wall kernel preserves the reference advancement, touching
check, and bisection for selector/top/bottom walls; the cam piston retains its
Python trajectory query. With all other paths held constant, 200 triangles
took 1.02 wall seconds with Python wall CCD and 0.81 with the compiled wall
kernel for 0.02 physical seconds; 500 triangles took 1.10 and 0.99 seconds
for 0.005 physical seconds. Batching circumcircle wall bounds and compiling
near-wall support gaps reduced the corresponding 200-triangle run from 0.78
to 0.70 seconds and the 500-triangle run from 0.92 to 0.88 seconds. The
small run-to-run differences between paired comparisons reflect timing noise.
Event counts, penetration, and energy residuals matched within each pair.
The full suite passes 88 tests; longer backend parity remains open.

The next search slice added an x-sorted swept circumcircle box search, with
identical conservative boxes in Python and Numba and deterministic body-index
ordering. On the same seed and machine, 200 triangles over 0.02 physical
seconds took 0.75 wall seconds with the grid and 0.55 with the compiled sweep;
500 triangles over 0.005 physical seconds took 0.91 and 0.52 seconds.
Event counts (43 and 14), zero penetration, and energy residuals matched.
These are short directional probes. The 92-test suite passes, including
randomized completeness and grid/sweep Carnot chronology checks at 200
triangles. Long-run
parity and interactive playback remain open.

The coarse batch pair-query slice calls the existing compiled polygon CCD
loop over sorted candidate pairs within one Numba invocation. It preserves
the scalar result format and has a `--pair-kernel scalar` switch for reference.
With the swept-box search fixed, 200 triangles over 0.02 physical seconds
took 0.53 wall seconds with scalar pair calls and 0.49 with the batch;
500 triangles over 0.005 physical seconds took 0.53 and 0.51 seconds.
Event counts (43 and 14), zero penetration, and energy residuals matched.
These short runs do not establish interactive playback or long-run parity.
The Carnot triangle default is now 96 particles, double the prior GUI count;
the configured radius rule keeps occupied area fraction fixed.

The larger-disc profile at 500 particles over 0.01 physical seconds assigned
about 1.77 of 2.89 profiled seconds to 156,000 wall query calls and 0.88 to
about 77,000 disc pair calls. A compiled conservative wall reach mask reduced
wall calls to about 6,300 in a comparable profile. With the same Numba
backend and scalar disc pair path, 200 discs over 0.1 physical seconds fell
from 1.96 to 1.47 wall seconds and 500 discs over 0.05 physical seconds
from 10.46 to 8.26. A compiled batch disc pair kernel then reduced short
200-disc runs from 1.53 to 1.08 seconds and 500-disc runs from 8.07 to
4.94 seconds. Different measurement pairs show run-to-run timing noise.
The disc area fraction was 7.54% of minimum cylinder area. The batch and
scalar runs had identical event kinds/participants (190 and 337 events),
maximum event-time differences below `3e-13`, position differences below
`5e-13`, velocity differences below `2e-11`, and heat/work differences at
roundoff scale. Both retained negligible penetration and first-law residuals.
The current suite passes 100 tests. These are short probes; long-run parity,
stationarity, and interactive playback remain open.

The controlled-shaft cam piston slice compiles the existing conservative
advancement and impact-time bisection for disc contacts in one batch. The
free-shaft trajectory stays on the Python path; `--cam-kernel python` selects
the controlled reference. With all other paths held constant, 200 discs over
0.1 physical seconds took 1.17 seconds with Python cam queries and 0.73 with
the compiled batch; 500 discs over 0.05 physical seconds took 5.28 and 2.91.
Event counts (190 and 337), penetration, and first-law residuals matched.
A 96-disc, 0.5-second run with 19 piston contacts had exact matching event
times, final state, and ledgers. The 102-test suite passes. A short post-change
500-disc profile assigns about 0.45 of 1.01 profiled seconds to the remaining
fixed-wall disc calls, 0.17 to penetration, and 0.12 to candidate search.
The replay option and long-run scientific gates remain open.

The fixed-wall disc slice batches conservative advancement and bisection for
the three constant-velocity infinite-line walls. `--wall-kernel python` keeps
the scalar reference; the free-shaft and controlled-shaft piston paths are
unchanged. With the same seed, 200 discs over 0.1 physical seconds took 0.70
wall seconds with Python fixed-wall queries and 0.30 with the compiled batch;
500 discs over 0.05 physical seconds took 2.94 and 1.58 wall seconds. The
500-disc pair had identical 337-event chronology, event times, final state,
heat, and work. The 96-disc 0.5-second trajectory also matched the reference.
The tests include moving, fast, touching, and separating direct disc-wall
queries. These short results still miss interactive playback. Penetration and
candidate construction are the next measured solver costs; compressed replay
and long-run scientific validation remain open. Finite membrane segments stay
on their separate segment solver; the full 105-test suite passes. A warmed
500-disc, 0.01-physical-second profile now assigns about 0.20 of 0.58
profiled seconds to penetration, 0.12 to swept-pair construction, and 0.08 to
decoding the disc pair batch. These are the next candidates for measured
optimization.

The disc-only penetration slice moves the ordered pair and infinite-wall
overlap check into one Numba call. Worlds with finite membrane segments or
polygon bodies retain the existing check, and `--penetration-kernel python`
keeps a reference path. On seed 123 with all other kernels fixed, 200 discs
over 0.1 physical seconds took 0.69 wall seconds with Python penetration and
0.23 with the compiled pass; 500 discs over 0.05 physical seconds took 3.22
and 1.24. The 500-disc run had identical 337-event chronology, event times,
final state, heat, and work. A new short 500-disc profile places about 0.12
of 0.40 profiled seconds in swept-pair construction and 0.09 in decoding
batched disc pair results. Those are the next measured targets. The run still
misses interactive playback; long-run validation and compressed replay remain
open. The full 107-test suite passes.

The swept-pair array slice keeps sorted candidate arrays through collision and
penetration queries. The public list-returning search remains available as a
reference. On seed 123 with the other kernels fixed, 200 discs over 0.1
physical seconds took 0.23 wall seconds before this change and 0.19 after;
500 discs over 0.05 physical seconds took 1.24 and 0.81. A same-process
500-disc comparison measured 0.81 seconds for the array path and 1.22 for
the original list path, with identical 337-event chronology, event times,
final state, heat, and work. A warmed 500-disc, 0.01-physical-second profile
now assigns about 0.027 of 0.264 profiled seconds to candidate construction
and 0.072 to decoding batched disc pair results. Pair-result decoding is the
next measured target. Interactive playback, long-run validation, and
compressed replay remain open. The full 108-test suite passes.

The deferred disc-pair decode slice keeps the raw batch result until wall and
portal queries finish, checks coincident-center failures in candidate order,
and constructs contacts only within the simultaneous window of the earliest
event. Pair-first order is retained when a wall contact has the same time.
On seed 123, the 200-disc, 0.1-physical-second run fell from 0.19 to 0.15
wall seconds, and the 500-disc, 0.05-physical-second run fell from 0.81 to
0.48. The 500-disc scalar pair reference had the same 337 event kinds and
participants; event-time differences were below `3e-13`, position below
`5e-13`, velocity below `2e-11`, and heat/work at roundoff scale. A warmed
500-disc, 0.01-physical-second profile now takes about 0.192 profiled seconds;
state deep copying accounts for 0.041, energy evaluation 0.034, and swept-pair
construction 0.027. Scheduler state copying and energy evaluation are the
next measured targets. These short runs remain below interactive playback;
long-run validation and compressed replay remain open. The full 110-test
suite passes.

The checkpoint and energy slice skips membership accounting in worlds without
portals, shares frozen cam geometry across checkpoint copies, and copies each
compensated ledger counter by value. Shaft, body, and mutable ledger states
remain independent across restore. On seed 123 with the Numba sweep backend,
200 discs over 0.1 physical seconds fell from 0.15 to 0.12 wall seconds;
500 discs over 0.05 physical seconds fell from 0.48 to 0.40. Event counts,
maximum penetration, and first-law residuals matched the prior short probes.
A warmed 500-disc, 0.01-physical-second profile assigns about 0.028 of 0.139
profiled seconds to candidate construction, 0.017 to penetration, and 0.014
to deep copying. The next session should reassess 200/500-triangle throughput
and longer-run parity before further disc-only tuning. Interactive playback,
scientific validation, and compressed replay remain open. The full 112-test
suite passes.

The triangle reassessment measured 200 particles over 0.02 physical seconds
at 0.35 wall seconds and 500 over 0.005 at 0.33 on the i7-11850H before the
new triangle diagnostic batch. The 200-triangle profile placed about 0.26 of
0.51 profiled seconds in penetration, mainly Python SAT separation. A batched
Numba current-state SAT gap path preserves body-major diagnostic order and
retains `--penetration-kernel python` as reference. On seed 123 with other
paths fixed, 200 triangles fell from 0.43 to 0.22 wall seconds over 0.02
physical seconds; 500 triangles fell from 0.41 to 0.31 over 0.005. Two
200-triangle seeds over 0.1 physical seconds and a 500-triangle seed over
0.02 had exact matching event chronology, times, final state, heat, work,
penetration, and CCD failures against the Python penetration path. A
200-triangle sampled checkpoint/restart over 0.1 also matched the direct run
exactly. A 200-triangle, 0.1-second comparison with scalar pair, Python wall,
and Python penetration paths matched all 226 event kinds/participants;
time differences stayed below `1e-15`, position below `4e-15`, and velocity
below `1.2e-12`. A post-change 200-triangle profile places about 0.11 of
0.31 profiled seconds in polygon wall queries, 0.08 in the polygon pair batch,
and 0.06 in penetration. These short probes remain far below interactive
playback. Multi-cycle backend parity, thermodynamic acceptance, and compressed
replay remain open. The full 115-test suite passes.

The controlled-shaft cam triangle slice batches the existing polygon-wall
advancement, touching checks, and bisection using the same quintic piston
profile. The free shaft retains its general solver, and `--cam-kernel python`
keeps the controlled reference. On seed 123 with other paths fixed, 200
triangles over 0.02 physical seconds took 0.23 wall seconds with Python cam
queries and 0.19 with the batch; 500 triangles over 0.005 took 0.31 and 0.22.
The 96-triangle, 0.5-second parity run included 24 piston contacts. Two
200-triangle seeds over 0.1 physical seconds and a 500-triangle seed over
0.02 matched the Python cam path exactly in event chronology, times, final
state, heat, work, penetration, and CCD failures. A post-change 200-triangle
profile assigns about 0.080 of 0.216 profiled seconds to the polygon pair
batch and 0.053 to penetration; polygon wall queries fell to about 0.018.
Polygon pair batching and diagnostic overhead are the next measured targets.
Interactive playback, multi-cycle validation, thermodynamic acceptance, and
compressed replay remain open. The full 117-test suite passes.

The swept-circumcircle slice computes the closest approach of each pair of
linearly moving centers over the search horizon and keeps pairs within the
sum of their fixed circumradii plus a geometry margin. This remains
conservative for rotating polygons. The compiled pair path uses this filter;
`--pair-kernel scalar` keeps the previous reach rule. At the initial seed-123
state and 0.05 horizon, the earlier filter retained 236 of 236 swept-box
candidates at 200 triangles and 1,061 of 1,062 at 500; the new filter retained
72 and 237. Short solver probes fell from about 0.19 to 0.14 wall seconds for
200 triangles over 0.02 physical seconds, and from 0.22 to 0.16 for 500 over
0.005. A same-process comparison holding the polygon CCD batch fixed measured
0.158 versus 0.171 seconds at 200, and 0.171 versus 0.227 at 500. Two
200-triangle seeds over 0.1 physical seconds and one 500-triangle seed over
0.02 matched the previous filter exactly in event chronology, times, state,
heat, work, penetration, and CCD failures. Randomized sampled trajectories,
tangency, and high-speed cases passed the conservative-mask check. Polygon
CCD and penetration diagnostics remain the largest measured costs. Interactive
playback, multi-cycle validation, thermodynamic acceptance, and compressed
replay remain open. The full 119-test suite passes.

The polygon CCD slice retains the full-witness scalar pair kernel as a
reference. In the compiled batch, conservative advancement and bisection use
the separating-axis gap; full witness geometry is evaluated at contact to
recover the normal and feature IDs. The scheduler also delays Python contact
construction until wall and portal queries establish the next event. Direct
SAT-gap and full-witness checks agreed within `2e-14` over randomized moving
triangle states. In same-process warmed comparisons holding the batched pair
search fixed, 200 triangles over 0.02 physical seconds took 0.116 wall
seconds with SAT-first CCD versus 0.134 with full witnesses; 500 triangles
over 0.005 took 0.127 versus 0.157. Two 200-triangle seeds over 0.1 physical
seconds and one 500-triangle seed over 0.02 matched exactly in event
chronology, times, state, heat, work, penetration, and CCD failures. The
benchmark now warms one complete event search to keep lazy compilation out
of solver timing. Its seed-123 short probes report 0.123 seconds at 200 and
0.135 at 500 triangles after warmup; compilation still contributes to
startup time. A warmed 200-triangle profile assigns about 0.054 of 0.168
profiled seconds to penetration and 0.045 to polygon pair CCD. Penetration,
polygon array packing, and multi-cycle backend parity are the next measured
targets. Interactive playback, thermodynamic acceptance, and compressed
replay remain open. The full 121-test suite passes.

A `cProfile` sample at 200 triangles before the fused transform found about
21,000 separation calls and 51,000 polygon transforms within a 0.005-second
physical interval. After compiling that full loop, a comparable profile put
polygon wall queries at about 0.10 of 0.21 profiled seconds and penetration at
about 0.07 seconds. Bounded wall/penetration checks reduced both; see the
measurements above. Further optimization remains subject to reference parity
and long-run collision gates. After the latest kernels, a 200-triangle
0.005-physical-second profile assigns about 0.018 of 0.064 profiled seconds
to penetration, 0.017 to pair calls, 0.016 to wall dispatch, and 0.010 to
grid construction. The cam piston still accounts for near-wall Python work.

After the fixed-orientation swept-axis CCD rejection, the matched isolated
three-seed, 0.1-physical-second triangle matrix at 200 and 500 particles is
saved in `runs/performance_matrix_swept_axis_2026-10-03.json`. All six event
counts, CCD-failure counts, maximum penetration values, and first-law
residuals exactly match the earlier `pair_scratch_after` matrix. Median warm
physical/wall throughput changed from 0.331 to 0.320 at 200 and from 0.0300
to 0.0297 at 500. These small short-window differences do not establish a
representative speed change; live interactive throughput remains unmet.

## 2. Establish thermodynamic evidence

- Extend the 18-cycle studies in
  `docs/CARNOT_THERMODYNAMIC_EVIDENCE_2026-09-30.md` to longer observation
  windows and the 96-triangle default after identifying a viable protocol.
  The original 40 runs had inconclusive gas-energy and hot-heat drift screens.
  The new three-seed, 36-cycle slow triangle jacket cohort has positive hot
  input in all four blocks for all seeds, but hot-heat drift remains
  inconclusive for all three and gas-energy shift is bounded for only one.
  The jacket raises slow-speed heat-gate readiness for both gas types.
- Investigate why gas enters hot branches above the hot reservoir temperature
  and why hot-wall heat has unstable signs. Separate cam timing, reservoir
  coupling, finite-size effects, and incomplete equilibration using branch
  and contact diagnostics. An optional cold-sector jacket now tests stronger
  reservoir coupling with the same predeclared heat and drift gates. Its
  matched ensembles improve heat-gate readiness from 0/3 to 2/3 for slow
  32-triangle runs and from 2/3 to 3/3 for slow 32-disc runs, but all drift
  screens in the 18-cycle cohort remain inconclusive. The 36-cycle extension
  alone does not resolve hot-heat drift. Next separate reservoir coupling and
  compression timing with matched protocol changes rather than repeating the
  same configuration or changing the gate to fit the data.
- With gate-ready ensembles, compare net measured efficiency and pressure-area
  work against the analytical Carnot reference as physical shaft speed falls.
  Report seed-level uncertainty and storage changes; do not infer convergence
  from a single favorable seed or a near-zero heat denominator.
- Extend pressure-area versus event-work comparison to free-shaft trajectories.
  Check heat/work signs and first-law closure for complete forward and reversed
  cycles over multiple seeds.

Gate: a reproducible multi-seed study supports or falsifies a slow-cycle trend,
with explicit stationary-behavior limitations, uncertainty, and no temperature
clamp. If the expected Carnot behavior does not emerge, document the physical
or numerical cause instead of publishing a measured efficiency claim.

The original ensemble documents a limited hot-heat speed trend. The optional
cold-jacket comparison improves heat-gate readiness for slow discs, while
stationarity remains inconclusive. This supports inadequate cold-sector
coupling as a contributor, but does not isolate it from cam timing or establish
stationary behavior. Section 2 remains open; the detailed evidence and
limitations are in `docs/CARNOT_THERMODYNAMIC_EVIDENCE_2026-09-30.md`.
The 36-cycle slow triangle jacket extension makes all three seeds heat-gate
ready, but all three hot-heat drift screens remain inconclusive. Its first 18
cycles exactly reproduce the earlier matched runs. Further unchanged long
runs are lower priority than a controlled physical change. A matched
16-disc hot-jacket test increased hot contact area roughly fourfold, but the
three paired hot-heat changes had mixed signs and all drift screens remained
inconclusive. An unequal-sector cam test also changed hot heat in mixed
directions without resolving drift. Section 2 remains open; further protocol
choices should use these matched results rather than repeat either setting.
The longer cold-sector test also reduced hot heat in all three matched seeds.
The existing matched 16- and 32-disc artifacts complete the particle-count
comparison; they do not satisfy the stationarity gate. Further unchanged
variants are lower priority than completing the open numerical, replay,
desktop, and release deliverables.
The 36-cycle, three-seed slow 32-disc cold-jacket extension reproduced the
first 18 cycles exactly and kept all three seeds heat-gate ready. All six
hot-heat and gas-energy drift screens remained inconclusive, so the stationary
positive-input and efficiency-trend checkpoint remains open. Details are in
the thermodynamic evidence document and the 2026-10-03 run artifact.
A matched 18-cycle 32-disc cold-jacket cohort at shaft speed 0.15 found higher
hot heat at speed 0.075 in all three seeds, but no faster seed passed the
positive-hot-heat block gate and all faster drift screens were inconclusive.
This supports a hot-input speed effect for this protocol without satisfying
the stationary-efficiency comparison gate.
The subsequent matched 32-disc both-jacket cohorts passed the positive-hot-
heat block gate at both speeds, while all drift screens remained inconclusive.
Faster hot heat exceeded slower hot heat in two of three paired seeds; see
`docs/CARNOT_GAP_FOLLOWUP_2026-10-04.md`.
Their 36-cycle extension again left all twelve drift screens inconclusive.
The first 18 cycles reproduced exactly, and the descriptive speed ordering
remained mixed by seed. Stop unchanged extensions and investigate the physical
protocol and screen before another long run; see
`docs/CARNOT_GAP_FOLLOWUP_36C_2026-10-04.md`.
The saved branch records localize the mean hot-entry excess to the preceding
compression: cold exit is 0.703–0.785, compression raises translational
temperature by 0.924–0.995, and mean hot entry is 1.636–1.766 against a 1.5
bath. A matched 18-cycle, 32-disc, both-jacket compression-timing pilot and
its decision rule are predeclared in `docs/CARNOT_NEXT_PROTOCOL_2026-10-04.md`.
That pilot failed its predeclared hot-entry reduction threshold at both shaft
speeds; all twelve drift screens remained inconclusive. No 36-cycle extension
of this variant is planned.
The next matched 64-disc, both-jacket pilot tests whether finite-count noise
explains the inconclusive drift screens. Its 18-cycle decision rule is fixed in
`docs/CARNOT_PARTICLE_COUNT_PILOT_2026-10-05.md`. The completed pilot found
zero bounded drift screens out of twelve, so the predeclared 36-cycle threshold
failed. Hot entry remained above the hot reservoir mean at both speeds. The
64-disc cohort does not establish stationary operation or justify an unchanged
extension. Section 2 stays open for physical model/protocol diagnosis and a
new controlled test; no measured efficiency convergence should be claimed.
The next matched dilution pilot varies only particle radius at fixed count,
with scale 0.5 and 0.25 at both speeds. Its decision rule is predeclared in
`docs/CARNOT_DILUTION_PROTOCOL_2026-10-05.md`. The completed pilot reduced
hot-entry temperature at slow speed but showed only a 0.015 mean reduction at
fast speed for quarter-scale radius, below the fixed 0.10 threshold. At fast
speed, lower compression heating was partly offset by warmer cold exit and
fewer hot contacts. No longer dilution run is planned. Further scientific
work needs a new physical protocol or a documented limitation, not another
unchanged cohort.

## 3. Finish numerical acceptance

- Expand triangle CCD coverage for grazing, nearly parallel, high-speed,
  rotation-only, and changing closest-feature contacts. Exercise hard
  near-simultaneous and multi-feature clusters.
- Test refinement and rollback on ambiguous contact cases; verify exact
  checkpoint/resume and isolated time reversal. Preserve diagnostic context
  for all unresolved cases.
- Run long dilute disc and triangle simulations over multiple seeds and sample
  cadences. Record penetration, CCD failures, first-law residuals, and event
  chronology. Keep `carnot_triangles` marked experimental until these pass.

Near-simultaneous three-disc and three-triangle regressions now place two
impacts half a time tolerance apart and check one elastic cluster,
sampling-cadence invariance, checkpoint replay, penetration, and energy
closure. The triangle case also compares Python and Numba event chronology
and final state. These pass; broader moving-boundary and multi-feature
coverage remains open.
Four targeted triangle contacts after cam-sector boundaries now agree among
the compiled cam batch, bounded scalar query, and unpruned scalar query.
A rotating triangle against a finite segment also retains its collision after
the nearest vertex changes during approach; direct, sampled, and checkpoint
replay histories match exactly. Longer moving-boundary trajectories and harder
multi-feature clusters remain open.
A symmetric two-triangle approach exposed a real second vertex-edge contact
at the first impact: resolving only the single CCD witness made that contact
penetrate immediately and exhausted horizon refinement. The impact resolver
now gathers approaching vertex-edge contacts within geometric tolerance and
passes the pair to the existing elastic cluster solver. The scene completes
two two-contact impacts, conserves energy, stays within penetration tolerance,
and returns to its starting state after velocity reversal. Direct, sampled,
and checkpoint-restored Python histories match; Numba agrees at numerical
tolerance. Broader multi-feature geometries and rollback refinement remain
open.
Two 32-triangle seeds now complete two controlled-speed cycles while bounded
and unpruned wall searches agree on same-state next-contact queries at all
eight branch midpoints per seed. Both bounded trajectories have zero CCD
failures, no material penetration, and small first-law residuals. For seed
123, a direct run, branch-midpoint sampled run, and replay from a second-cycle
checkpoint exactly match event history, final state, and ledger. Independent
bounded and unpruned trajectories differ at roundoff-level event times early
and later separate in event chronology, so exact long-run independent-history
parity is not claimed. Harder cluster and broader seed/cadence checks remain
open.
Fixed-orientation polygon pair CCD now rejects a query only when swept
projections remain separated on one polygon edge axis for its full horizon.
A near-parallel triangle miss that previously caused 3,109 scheduler horizon
refinements now returns no collision with zero refinements; a nearby grazing
impact retains its original time in Python and Numba, including at high
speed. Rotating pairs retain conservative advancement. Broader grazing and
changing-feature acceptance remains open.
Two additional 96-triangle seeds completed two cycles with exact sampled and
checkpoint replay at distinct cadences and matching branch-midpoint same-state
backend queries; see `docs/CARNOT_GAP_FOLLOWUP_2026-10-04.md`.
The next batched acceptance probe covers dilute 32-disc and 32-triangle box
trajectories for 50 physical seconds at seeds 128–129. It compares the direct
trajectory with two sample cadences and a midpoint checkpoint replay, and
records event chronology, penetration, CCD failures, and first-law residual.
This extends the long-run sampling evidence; the harder moving-boundary and
multi-feature cases above remain separate requirements.
The completed batch passed all four trajectories, all 156 unit tests, and
scientific validation. Each trajectory recorded 4,657–5,562 events with
monotone chronology, zero CCD failures, and absolute first-law residual below
`4e-14`. Penetration was at most `1.6e-14`, versus `1e-10` geometry tolerance.
Both 0.37 s and 1.11 s sampling schedules and midpoint checkpoint replay
matched the direct event history, ledger, and final body state exactly.
The compact report is
`runs/gap_followup_2026-10-04/numerical_dilute_50s_2026-10-05.json`.
The report's sampled CCD refinement totals include the additional replay
work; they are computational counters, not a sampling-dependent event change.
The remaining targeted numerical batch runs three independent checks in
parallel: two fast-shaft cycles at 32 triangles across seeds 130–132 with
moving-wall Python/Numba next-contact comparisons; 300 separated rotating
triangle-pair CCD queries checked against both backends and dense overlap
sampling; and seven offsets around a symmetric two-feature triangle impact
with sampling, checkpoint replay, and backend comparison. Any mismatch will
be investigated before closing section 3.
The first hard-contact batch completed the moving-boundary and cluster jobs:
all three seeds finished two fast-shaft cycles with zero CCD failures,
penetration, or branch-midpoint Python/Numba query disagreements; seed 130's
sampled and checkpoint-restored trajectories matched direct. All seven
two-feature offsets matched direct, sampled, checkpoint, and compiled results,
with zero CCD failures or penetration. The random pair job completed its
queries but failed while serializing a NumPy integer, before reporting a
scientific verdict. The successful reports are saved as
`runs/gap_followup_2026-10-04/numerical_hard_moving_boundary_2026-10-05.json`
and `numerical_hard_offset_clusters_2026-10-05.json`. The next launcher retries
only the random pair report and adds an independent tolerance-convergence
probe; it does not repeat the passed jobs.
The retry passed both probes. All 300 separated rotating triangle pairs had
matching Python/Numba status and contact time within `1e-8`, no observed
overlap classified as no collision, and no indeterminate result. Rotation-only
contact, a nearly parallel hit, and a nearby miss kept their classifications
at geometry tolerances `1e-8`, `1e-10`, and `1e-12`; both backends agreed and
the hit times converged. Reports are
`runs/gap_followup_2026-10-04/numerical_hard_random_pairs_2026-10-05.json`
and `numerical_hard_tolerance_2026-10-05.json`. Together with the 50-second
dilute cohorts, earlier high-count triangle cohorts, and strict diagnostic
handling for unsupported simultaneous policies, these satisfy section 3's
declared acceptance set. The numerical gate is closed for that set; this is
not a proof for arbitrary geometries or the brief's aspirational million-impact
engineering target, which belongs to release evidence.

Gate: no unreported overlap or skipped event, and event chronology is
sampling-independent within declared tolerances across the acceptance set.

## 4. Complete desktop controls and presentation

- Add reset-to-seed, fixed-duration step, one-branch step, playback-speed and
  interactive physical shaft-speed controls, preset selection, checkpoint,
  export, and a clear healthy/refining/paused state.
- Add camera zoom, pan, and selection; an expandable energy/ledger plot; plot
  export, saved configurations, recording, and explanations/inspectors needed
  to explore the apparatus.
- Check real-display responsiveness and minimum-window layout during slow
  physics computations and dropped render frames.

Gate: GUI and headless runs with the same seed and command schedule agree, and
all controls remain responsive without changing physical event horizons.

The first control pass adds fixed-duration and next-branch stepping, seed
reset, preset switching, target playback rate, and a separate physical
shaft-speed control. The latter logs an explicit motor-work intervention and
updates the controlled phase origin; checkpoint restore preserves that origin.
The apparatus view now has zoom, pan, selection, and fit-to-scene. Focused
offscreen UI and physics tests pass.
The real xcb desktop run matched the headless position, velocity, event count,
and energy residual after the same shaft-speed and branch-step schedule. The
1.136 s branch step allowed 43 GUI timer callbacks; pan, selection, preset
switch, seed reset, and the 850×520 window also passed. All 153 unit tests and
scientific validation passed. The compact result and screenshot are in
`runs/gap_followup_2026-10-04/desktop_controls_acceptance.json` and
`desktop_controls.png` in the same folder. The 850 px screenshot showed toolbar
overflow, so the toolbar is now split into two rows. Configuration and checkpoint
save/restore, run export, and three PNG plot exports pass a focused offscreen
round trip. An expandable energy/ledger plot is also present. Real-display
verification of these new actions and layout, recording, explanatory inspectors,
and dropped-frame responsiveness are still required before this gate closes.
The next desktop batch also tests timestamped screen-frame recording.
Its first attempt passed 154 unit tests and scientific validation. Desktop
checkpoint restore, run export, three plot exports, toolbar fit, and ten
recorded frames completed, but a GUI timer-count assertion failed before the
script printed its timing counts. Recording timestamps showed irregular frame
spacing up to about 0.63 s while the branch step ran. PNG compression now runs
in a bounded background writer, and the repeat batch will report GUI timer
gaps, capture time, and dropped recording frames before judging responsiveness.
The attempt's compact diagnosis and final screenshot are saved as
`runs/gap_followup_2026-10-04/desktop_files_attempt1.json` and
`desktop_files_attempt1.png`.
The second real-display attempt confirmed the two-row toolbar fits 850×520,
checkpoint restore and exports work, and the full 154-test and scientific
suites pass. The isolated standard benchmark also passed. Recording still
blocked GUI timer callbacks: nine frames took up to 0.418 s each to capture,
leaving only ten 20 ms timer callbacks over a 2.386 s branch step and a maximum
callback gap of 0.436 s. The asynchronous PNG writer therefore did not solve
the expensive QWidget capture. The second attempt is preserved as
`runs/gap_followup_2026-10-04/desktop_files_attempt2.json` and its screenshot.
The third batch compared physics steps without and with recording through an
X11 screen capture path, while keeping the paused/refining state visible after
file actions.
The third xcb batch passed all five jobs: 154 tests, scientific validation,
environment capture, the isolated standard benchmark, and desktop acceptance.
At 850×520 there was no toolbar overflow. The same checkpoint restored twice;
run and three plot exports succeeded. The 10.472-physical-second branch step
took 2.438 s without recording and 2.484 s with recording; corresponding
20 ms GUI timer maxima were 46.3 ms and 51.8 ms. The X11 capture path saved
26 frames with no drops, a maximum capture cost of 18.7 ms, and a final frame
at the branch endpoint. The persistent state returned to PAUSED. The compact
result, screenshot, and environment are in
`runs/gap_followup_2026-10-04/desktop_files_acceptance.json` and
`desktop_files_acceptance.png`. A focused offscreen test additionally verifies
that dropping intermediate render frames leaves the final gas-box trajectory
and event count identical to headless, while click inspection displays live
particle and apparatus details. This closes the desktop gate for the measured
controlled-shaft desktop scope; other hardware and display backends have not
been measured.

## 5. Release acceptance

- Run a clean Ubuntu installation and the full scientific and GUI suite.
- Publish long-run conservation, multi-seed thermodynamic, collision, 200/500
  particle performance, and screenshot artifacts with hardware and backend
  details.
- Map every Carnot requirement in the original implementation brief to a
  test, benchmark, screenshot, or exported run artifact.

Gate: all applicable requirements have reviewable evidence. The experimental
label remains until the numerical and scientific gates above pass.
The 36-cycle study batch also passed 144 unit tests and all scientific
validation checks; GUI and release acceptance remain open.
The 2026-10-05 fresh-venv run on Ubuntu 24.04.5 LTS and Python 3.12.3
installed the GUI and Numba extras, passed 156 unit tests and all 16 scientific
checks, and confirmed that a headless import did not load Qt. This covers the
local fresh-environment part of installation acceptance; a clean OS image and
the complete requirement-to-evidence map remain open.
An initial Carnot-specific map of the original brief and its remaining gaps is
in `docs/CARNOT_RELEASE_EVIDENCE_MAP_2026-10-05.md`. Other experiments in the
brief still need release mapping.
The prepared clean-image check is described in
`docs/CARNOT_CLEAN_UBUNTU_PROTOCOL_2026-10-05.md`. It has not been run; this
host has no Docker or Podman executable, so evidence must come from a separate
clean Ubuntu system.

## 6. Add precalculated, smooth replay

User-requested playback option: calculate a run ahead of time, save it in a
compact, readily compressible format, then replay it in real time without
waiting for collision solving. This is a separate delivery path from live
simulation; continue the current performance, numerical, and scientific
gates above rather than substituting replay for them.

- Specify a versioned run format with seed/config/backend metadata, geometry,
  exact event and ledger records, and timestamped particle/apparatus states.
  Evaluate chunked columnar arrays with compression, periodic keyframes, and
  bounded-error deltas; quantify file size and maximum reconstruction error.
- Add an offline precalculation command that streams chunks to disk so long
  runs do not need all frames in RAM. Store enough state for seeking and
  responsive playback at user-selectable speed.
- Build replay in the desktop view with smooth time-based interpolation of
  positions, angles, and apparatus motion. Show statistics from recorded
  values at the correct timestamps, with event-aware handling of discontinuous
  velocities and branch changes.
- Benchmark 200/500-particle replay for sustained frame rate, seek latency,
  load time, memory, and compressed size on declared hardware. Verify saved
  samples and observables against the source run, including after a seek.

The v1 archive now streams float64 frame chunks to compressed columnar files,
records exact events and final ledger state, and supports random frame access
and seeking. Saved-frame reconstruction is exact. A desktop replay window
provides time-based playback, speed, seek, and timestamped recorded stats;
interpolation uses event times for branch and velocity changes but remains
approximate between frames with collisions. Short 200/500-disc size,
cold-seek, and offscreen render probes are in `docs/REPLAY_FORMAT.md`.
Long-run storage has now been measured; exact between-event reconstruction
and smooth real-display acceptance were subsequently addressed for measured
cases. The earlier archive writer retained the solver's full event history.
Ten-physical-second 200/500-disc archives now quantify longer-window storage
and offline write cost; see `docs/CARNOT_GAP_FOLLOWUP_2026-10-04.md`.
The real-display `xcb` check reached 60.7–61.1 painted frames/s over ten
physical seconds, with 95th-percentile intervals near 16.7 ms. Each archive
still had three gaps over 25 ms, reaching 149–183 ms, and cross-chunk seeks
took up to 55 ms. The rate target is met on this machine, but the smoothness
gate remains open; see `docs/CARNOT_GAP_REPLAY_DISPLAY_2026-10-04.md`.
The reader now warms the event index and first three chunks before window
display and prefetches later chunks with a bounded cache. Focused replay tests
and the 145-test full suite pass. A repeated real-display check on the same
ten-second archives reached 62.5 paints/s at both counts, with no interval
above 18 ms and maximum measured seek of 25.2 ms. This closes the observed
paint-gap problem for those archives at 1×. Longer random seeks and
between-frame error remain open; see
`docs/CARNOT_GAP_REPLAY_DISPLAY_2026-10-04.md`.
The one-second midpoint check found a maximum interpolated disc-position
error of 1.73 and 3.39 particle radii at 200 and 500 discs, despite exact
saved-frame agreement. Event-aware disc reconstruction is needed before the
replay-fidelity gate can close.
Event-aware disc reconstruction is now implemented with backward-compatible
reading of the existing archives. Independent 200/500-disc source comparisons
show roundoff-level position error and exact velocity at 60 midpoints over one
second. Both ten-second desktop runs still average about 62.5 paints/s, but
each had one 39–41 ms gap away from a chunk transition. The 146-test suite
passes. The current gate remains open for consistent smoothness, broader
fidelity, uncached seeking, and triangle angular replay; see
`docs/CARNOT_GAP_REPLAY_DISPLAY_2026-10-04.md`.
The next eight-job batch localized the 1× long paint gaps to Python garbage
collection, with two measured collection pauses of 24–29 ms. Repacked
19-chunk seek checks stayed below 10 ms uncached at both disc counts; 1× and
4× desktop runs all averaged about 62.4 paints/s. New 32/96-triangle
midpoint checks found errors up to 0.52/1.02 particle radii and π radians.
These measurements and remaining caveats are in
`docs/CARNOT_GAP_REPLAY_DISPLAY_2026-10-04.md`.
Playback now suspends cyclic garbage collection only while playing, and new
triangle archives record contact points for event-aware angular replay.
The acceptance batch ran 200/500-disc desktop playback at 1× and 4× for ten
physical seconds, plus 32/96-triangle playback at 1× for two physical
seconds. All six held about 62.5 paints/s without a gap above 25 ms, and
garbage collection was enabled again afterward. Independent source comparisons
at 120 midpoints for each triangle count found maximum position errors below
`7e-16`, angular errors below `7e-14` radians, and angular-velocity errors
below `1.2e-11`. All 147 unit tests and scientific validation passed. These
single-seed, short triangle checks do not yet establish longer replay behavior,
seek/statistic agreement after arbitrary seeks, or bounded precalculation
memory. See `docs/CARNOT_GAP_REPLAY_DISPLAY_2026-10-04.md`.
Seed-124 five-second 32/96-triangle comparisons now also pass at 300
midpoints each with roundoff-level motion error, exact saved scalar statistics
and ledger, and 62.5 paints/s on the desktop without a gap above 25 ms.
Random frame/ledger seeks at 200/500 discs matched the reference reader.
Archive-generation memory and broader high-count sustained playback remain
open; see the replay display report.
The replay writer now drains committed events to JSONL at every saved frame
and builds the original version-1 index using temporary disk-backed arrays.
The ten-second 200/500-disc regenerated archives match the old event logs by
SHA-256 and all saved frames and ledgers exactly. Peak generation RSS fell to
about 150/155 MiB. All four 1×/4× desktop runs of the genuine 19-chunk
archives reached about 62.5 paints/s without a gap above 25 ms. The 148-test
suite and scientific checks passed. See the replay display report. The
default controlled-cam point interpolation now has a conservative analytic
bound of `3.085e-6` for discs and `5.641e-6` for triangles at speed 0.15 and
60 saved frames/s; selector-shoe jumps use the exact logged branch time.
The 10.6-second 200/500-triangle archives now pass source and geometry
comparisons across that transition. The first 500-triangle 4× paint count
covered less wall time than expected, so five display-only timing runs measured
startup, ticks, paints, and play-to-end duration. Four 4× repeats each counted
166 paints, reached the end in about 2.657 wall seconds, and had no gap above
17.5 ms. The earlier low count did not recur and its cause is unknown.
The section 6 replay gate is accepted for the tested controlled-shaft Carnot
presets at 200/500 particles on the declared desktop, with storage, observed
fidelity, and analytic apparatus bounds documented in the replay reports.
Longer archives, free-shaft acceleration, and arbitrary operating-system
scheduling remain outside this measured acceptance set.

Gate: a precalculated run replays at the display target frame rate without
solver stalls, with documented storage/error bounds and faithful statistics.
