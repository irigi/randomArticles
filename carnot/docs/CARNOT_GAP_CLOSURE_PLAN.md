# Carnot gap closure: remaining work

Updated: 2026-09-29

The gap is not closed. The application has a working experimental apparatus,
plots, diagnostics, event accounting, and a reproducible speed study. The full
change history and completed tests are in `docs/PROGRESS.md`; this file tracks
only work still to do. The current 16-disc, 18-cycle probes have inconclusive
stationarity screens and no gate-ready measured efficiency. The collision engine
also needs broader long-run validation before `carnot_triangles` can lose its
experimental label.

## 1. Make larger particle runs practical

Selected route: **A+B**, conservative spatial search plus optional Numba
kernels. The implemented stages keep a fixed,
dilute occupied area fraction as count changes (triangles 3.12%, discs 7.54%
of the minimum cylinder area); a swept circumcircle grid replaces the dense
pair search from 32 particles; the numerical reach filter and full polygon
pair conservative advancement have Numba paths; and cam branch boundaries
are cached. Swept circumcircle bounds prune wall queries, and Numba kernels
now handle fixed-velocity polygon wall CCD and batched wall penetration
bounds. The unpruned wall and all-pairs/Python paths remain selectable for
reference comparisons. This is a measured partial speedup,
not completion of the 200/500-particle performance gate.

- Extend the reproducible Carnot benchmark for discs and triangles. Record warm
  wall time per physical second, events per wall second, memory, startup cost,
  particle count, seed, shaft speed, CPU, and backend. Profile representative
  16, 32, 64, 128, 200, and, when physically feasible, 500-particle cases.
- Profile the remaining scheduler and geometry costs at 200/500 particles.
  The fixed-velocity wall loop and wall penetration bounds are compiled;
  the cam piston still uses Python motion queries. A short 200-triangle
  profile now spreads time across pair-call orchestration, swept-grid
  construction, wall dispatch, and penetration. Test coarse batch pair
  queries or a cheaper grid only when they preserve event order; profile
  larger-count discs before compiling their TOI. Keep `fastmath` disabled
  and the Python path.
- Improve particle placement at 200/500 particles if setup time warrants it.
  Confirm the radius rule, density, and finite-size effects in every report.
- Keep the swept-grid bounds conservative for rotating bodies; add targeted
  high-speed, grazing, boundary, and moving-piston parity cases. Extend
  bounded-wall parity across longer moving-piston trajectories and seeds.
- At every stage compare backend event chronology, contacts, heat/work ledger,
  penetration, checkpoint/restart, and failure diagnostics against the Python
  reference across fixed seeds, triangle edge cases, and multiple sample
  cadences. Measure speed only after the physics gate passes.

Gate: the chosen backend provides a measured improvement on representative
Carnot runs, with unchanged scientific results within declared tolerances. Aim
for 200 triangles at an interactive achieved playback rate, then validate 500
at a documented density and hardware configuration. Publish the benchmark and
any cases that miss the target. Compilation time and fallback behavior count
as part of the result.

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

## 2. Establish thermodynamic evidence

- Run larger independent seed ensembles with enough complete cycles for the
  predeclared transient and cycle-batch screens. Report gas-energy and hot-heat
  shifts, uncertainty, thermal contact counts, and numerical health per seed.
- Check whether hot heat is consistently positive in contiguous post-transient
  blocks. Revisit cycle speed, reservoir coupling, finite-size effects, and
  observation length using measured branch data if it is not.
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

## 5. Release acceptance

- Run a clean Ubuntu installation and the full scientific and GUI suite.
- Publish long-run conservation, multi-seed thermodynamic, collision, 200/500
  particle performance, and screenshot artifacts with hardware and backend
  details.
- Map every Carnot requirement in the original implementation brief to a
  test, benchmark, screenshot, or exported run artifact.

Gate: all applicable requirements have reviewable evidence. The experimental
label remains until the numerical and scientific gates above pass.
