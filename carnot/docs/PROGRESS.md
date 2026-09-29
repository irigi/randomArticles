# Progress record

## 2026-09-28 — reference implementation

Implemented the installable package, state arrays, deterministic event loop,
disc and rotating convex collision queries, finite material segments, elastic
and thermal policies, energy-step portals, measurements, exports, CLI, optional
Qt view, and validation suite.

Added these runnable experiments: thermal disc box, triangle equipartition,
controlled/free and forward/reversed Carnot apparatus, hard labyrinth,
energetic labyrinth, and selective membrane.

Validation evidence on CPython 3.12.3 and NumPy 2.4.6:

- 500 randomized isolated binary impacts: maximum normalized energy,
  momentum, and angular-momentum residuals below `2e-14`.
- Analytic disc/disc, tangent, and wall impact times within `1e-12`.
- Thermal flux exponential-mode mean within 0.004 of the target in 20,000
  samples (gate 0.025).
- Controlled Carnot short-run energy residual about `1e-14`; free-shaft run
  about `3e-15`.
- Hard and energetic labyrinth and selective-membrane smoke runs close their
  ledgers at floating-point scale.
- Ten deterministic unit tests pass, including checkpoint/resume equality.

The test suite is intentionally a reference-engine gate. Long-run statistical
acceptance, large particle benchmarks, mobile compound hosts, general contact
clusters, and complete GUI polish remain open and are listed in the README.

## Gap closure session 1 — triangle contacts and strict failure handling

- Replaced guessed polygon contact points with closest vertex/edge witness
  points and stable feature IDs. Added exact triangle support geometry for
  infinite-line cylinder and piston walls.
- Continued CCD search after a separating contact so a different rotating
  feature can collide later in the same interval. Increased the advancement
  budget for near-touching polygon and polygon-wall searches.
- Added initial-state and post-advance penetration checks. A severe overlap or
  indeterminate query pauses the run and retains a pre-failure checkpoint.
  The optional GUI writes `diagnostic.json` and `state.npz` to a temporary
  diagnostic directory on failure.
- Fixed invalid initial disc/host overlaps in the selective-membrane preset.
- Added regressions for the seed-123 Carnot triangle overlap, an exact
  triangle-wall time of impact, a high-speed wall impact, polygon witness
  locations, and a rotation-only impact whose interval endpoints are both
  disjoint. All 15 unit tests and the scientific validation suite pass.
- A 48-triangle controlled Carnot run to physical time `1.0` completed 338
  events with zero measured penetration and energy residual `3.24e-14` on the
  available Python 3.12/NumPy environment.

Remaining before milestone 2 acceptance: broad randomized and grazing CCD
cases, certified handling of near-simultaneous clusters, long-run
multiple-seed penetration testing, and sampling-cadence invariance. The GUI
apparatus and statistics panel remain part of later gap-closure phases.

## Gap closure session 2 — sampling-independent interval scheduling

- The scheduler now searches each internal interval from a fixed physical
  origin. Intermediate `advance_to()` calls materialize a snapshot from that
  origin and preserve the pending event, so render cadence does not alter the
  next collision query. Checkpoints retain this pending interval, and restore
  truncates the event history to the checkpoint's event count.
- Indeterminate CCD queries now enter a bounded 32-attempt horizon refinement
  loop. Failure diagnostics retain the query reasons and attempted horizons.
  The one-step command uses the same refinement path.
- Added tests for irregular snapshot cadence, checkpoint/resume during a
  pending interval, and repeated CCD refinement. All 18 unit tests pass.
- The seed-123 16-triangle reproduction to `t=0.2` now has exactly identical
  event times, positions, and velocities for one call and irregular calls.
  The 48-triangle run to `t=0.5` matches exactly for one call and 50 calls:
  161 events, zero measured penetration, and energy residual `-4.52e-14`.
  Four 16-triangle seeds reached `t=0.5` with no measured penetration.

Phase 2 remains open for near-simultaneous contact clusters, rollback with
geometry-accuracy refinement, and time-reversal coverage. Phase 1 still needs
polygon/finite-segment CCD and the wider isolated-impact test matrix.

## Gap closure session 3 — simultaneous stationary elastic contacts

- The scheduler now gathers contacts inside the earliest-time tolerance window
  and groups contacts that share a body. Independent groups remain independent.
- Small stationary elastic groups (up to six contacts) use a rank-aware active
  set solve. The solver accepts nonnegative impulses only when post-impact
  normal velocities are admissible and kinetic energy closes. It reports the
  matrix rank and residual in event metadata. Unsupported simultaneous thermal,
  portal, or moving-mechanism groups pause with a retained checkpoint.
- Checkpoints now retain event history for restoring a later checkpoint after
  branching from an earlier one. Cluster counters also restore correctly.
- Added exact simultaneous three-disc and corner-wall regressions, a duplicate
  wall rank-deficiency case, checkpoint branch restoration, and a thermal
  corner failure case. All 23 unit tests and the scientific suite pass.
- The 48-triangle Carnot run to `t=0.5` still has 161 identical events for one
  call and 50 sample calls, zero measured penetration, and zero clusters in
  that particular run.

Phase 2 still needs near-simultaneous time refinement, rollback with geometry
accuracy refinement, larger multi-seed runs, and time-reversal validation.
The present cluster solver covers stationary elastic contacts only.

## Gap closure session 4 — finite polygon-wall geometry

- Added exact closest-feature queries for a convex polygon against a finite
  segment, including endpoint, edge, crossing, and contained-segment cases.
- Triangle impacts on capsule-like finite walls now use conservative
  advancement with the actual polygon vertices and edges, not the triangle
  circumradius. The contact point comes from the polygon and capsule surface
  witnesses, allowing the correct off-center angular impulse.
- Initial-state and committed-state penetration checks now include finite
  segments for both discs and polygons.
- Added face-impact, clear-miss, rotation-only, high-speed, crossing, and
  initial-overlap regressions. All 28 unit tests and the scientific validation
  suite pass. Short seed-123 runs of both labyrinth presets, selective membrane,
  and Carnot triangles reach `t=0.5` without material penetration or energy
  residual above floating-point scale.

Phase 1 still needs broader randomized grazing and isolated triangle-impact
coverage. Phase 2 still needs time refinement for near-simultaneous contacts
and rollback with geometry-accuracy refinement.

## Gap closure session 5 — stable face normals and triangle reversal

- Vertex-edge and edge-vertex polygon contacts now take their normal from the
  contacted edge. The previous normal used the difference of two witness
  points only about `1e-10` apart, which made return impacts sensitive to tiny
  coordinate errors.
- Added a `reverse_particle_velocities` command for controlled time-reversal
  probes. In an oblique isolated two-triangle collision, the round-trip angle
  error fell from about `6e-7` to `6e-11` radians at default tolerances.
- The scientific suite now includes 200 seeded approaching triangle scenes.
  It found 198 impacts; the maximum normalized energy, linear-momentum, and
  angular-momentum residuals were `1.05e-15`, `4.45e-16`, and `1.12e-16`,
  with zero measured post-impact penetration.
- All 29 unit tests and the scientific suite pass. Four 16-triangle Carnot
  seeds reached `t=1.0` without measured penetration. The seed-123
  48-triangle run to `t=0.5` still had 161 identical events across one call
  and 50 sample calls.

The symmetric two-triangle scene with multiple successive feature contacts
still has a much larger reversal error and needs manifold/chronology work.
Broader grazing and nearly parallel cases, rollback refinement, and
near-simultaneous time refinement remain open.

## Gap closure session 6 — shared free-shaft trajectory

- Piston CCD, selector branch queries, committed shaft integration, and piston
  velocity now use one trajectory interface. Free-shaft prediction includes
  spring and load acceleration instead of extrapolating at constant speed.
- A zero-load shaft spring uses the exact harmonic trajectory, preserving its
  mechanical energy at each internal step. The load-bearing path retains its
  velocity-Verlet step; its predicted phase now matches the committed phase.
- Moving-piston CCD uses the trajectory's speed bound and a finer local root
  tolerance. Simulation updates the mechanism clock directly, removing the
  preset loader's advance-method wrapper.
- Added an accelerating-piston collision regression and an exact spring
  trajectory/energy regression. All 31 unit tests and the scientific suite
  pass. Seed-123 controlled and free Carnot triangle runs to `t=0.5` finish
  with zero measured penetration and floating-point-scale ledger residuals.

Phase 3 remains open for reflected piston inertia, exact branch-boundary
events, load-bearing conservative integration and full-cycle work checks.
The apparatus UI and instruments are still pending.

## Gap closure session 7 — exact cam branch events

- The scheduler now creates a physical `branch_transition` event at each cam
  sector boundary, before searching further contacts. The transition carries
  unwrapped phase, boundary index, and source/destination branch metadata.
- Controlled forward and reversed cycles both produce four transitions per
  revolution. The spring-driven free shaft also detects a return crossing
  after its turning point. Event times are independent of snapshot cadence.
- Selector queries retain the cam's right-continuous phase convention.
  Branch transitions split internal intervals so contacts after a boundary
  use the new branch state. The one-collision command skips branch records.
- All 33 unit tests and the scientific suite pass. Short controlled and free
  Carnot triangle runs retain zero measured penetration and floating-point
  ledger closure.

Phase 3 remains open for reflected piston inertia, load-bearing conservative
integration, and complete-cycle first-law checks. The apparatus UI and live
statistics are still pending.

## Gap closure session 8 — load-bearing shaft reversals

- The constant-inertia free shaft now follows exact piecewise motion under a
  linear spring and constant opposing load torque. It splits at angular-speed
  zeroes and sticks when spring torque cannot overcome the load.
- The load counter uses total angular travel, including reversals. Piston CCD,
  branch scheduling, and committed motion query the same trajectory.
- Added reversal, split-step, and simulation-ledger regressions. A seeded
  500-case shaft probe had worst absolute first-law residual below `8e-15`.
  All 35 unit tests and the scientific suite pass. Both controlled and free
  seed-123 Carnot triangle runs to `t=0.5` retain zero measured penetration
  and floating-point-scale energy residuals.

Reflected piston inertia, controlled-mode spring/motor accounting, and
complete-cycle first-law checks remain open in Phase 3. The apparatus UI and
live statistics are still pending.

## Gap closure session 9 — reflected piston inertia in isolated mechanisms

- Added the cam-linked piston mass to the shaft's effective inertia,
  `M(phi) = I + m_p f'(phi)^2`, and used it in kinetic energy, shaft speed,
  free piston recoil, and CCD speed bounds.
- Force-free variable-inertia motion now advances by the one-coordinate
  kinetic arc length and recomputes momentum at the new cam angle. A split
  step matches a direct step, and energy closes through a particle-piston
  collision with reflected mass.
- Controlled shafts now record smooth motor work from changing kinetic and
  spring energy, plus work at each piston impact. Load output is tracked
  separately. All 39 unit tests and the scientific suite pass.

Reflected piston mass is an optional shaft parameter. The loaded free Carnot
preset still uses zero piston mass because the coupled nonlinear spring/load
flow has not been implemented; constructing that unsupported combination now
raises a clear error. The next Phase 3 step is a conservative nonlinear
partial-step solver for that combination, followed by complete-cycle ledger
checks. The apparatus UI and live statistics are still pending.

## Gap closure session 10 — coupled nonlinear shaft flow

- Implemented an energy-conserving partial-step map for free shafts with
  reflected piston inertia, a linear spring, and opposing load torque. The
  scalar solve uses the energy identity and a symmetric discrete-gradient
  position equation. Steps split at shaft reversals and accumulate total
  angular travel for load work.
- Piston CCD, branch-boundary searches, and committed motion evaluate the
  same partial-step map. Both Carnot presets now include piston mass `0.5`.
- All 42 unit tests and the scientific suite pass. The coupled reversal and
  accelerating-piston collision tests close energy to floating-point
  precision. A 500-case randomized probe had
  worst absolute first-law residual `1.07e-14`. The 16-triangle free preset
  has identical events and state across direct and irregular sampling. At a
  `0.05` internal nonlinear step cap, 48-triangle controlled/free runs to
  `t=0.5` had zero measured penetration and residuals below `2e-14`.
- Reduced the nonlinear step cap to `0.01` after a trajectory-refinement
  comparison; all 42 tests still pass at this setting. The refined solver is
  slower and has not yet undergone a long complete-cycle acceptance run.

Phase 3 still needs documented trajectory convergence against an independent
integrator, complete-cycle first-law checks, and named apparatus components.
Phases 4 and 5 still need the visible apparatus, worker, plots, and instruments.

## Gap closure session 11 — independent trajectory and cycle checks

- Added a conventional RK4 reference for the canonical variable-inertia
  Hamilton equations with spring and opposing load. Over a smooth `0.08` time
  interval, the coupled solver's maximum phase/momentum error fell from
  `1.041e-3` at step `0.04` to `2.609e-4` at `0.02` and `6.526e-5` at `0.01`,
  consistent with second-order convergence. Two RK4 reference step sizes
  agreed within `1e-12`.
- Added complete forward and reversed controlled-cycle checks of branch
  transitions, motor work, load output, stored energy, and first-law residual.
  A complete free cycle with finite reflected piston mass also closes stored
  energy against load output. These checks use an isolated apparatus without
  gas collisions or thermal exchange; they do not establish thermodynamic
  cycle convergence.
- All 45 unit tests and the scientific validation suite pass. The free-cycle
  phase at a time found from a single uninterrupted trajectory differed from
  the interval-scheduled phase by `1.57e-7` radians, reflecting finite-step
  trajectory error. The cycle test states this tolerance explicitly.

Phase 3 still needs named geometry/energy state for the piston, flywheel,
spring, cam follower, selector, and load, plus full-cycle gas/thermal first-law
checks. Phases 4 and 5 still need the visible apparatus and instruments.

## Gap closure session 12 — named apparatus snapshot state

- Added immutable component records for the cylinder, piston, piston cam,
  follower, selector cam and shoe, thermal wall, flywheel and spoke, torsion
  spring, and load. Simulation snapshots expose their scene-coordinate paths,
  motion, branch state, energy, and cumulative load output.
- The piston path and speed use the exact cam functions used by moving-wall
  collision detection. The piston and flywheel kinetic energies plus spring
  potential energy reproduce the mechanism's stored energy without assigning
  extra energy to schematic cams or followers. Both cam track paths rotate
  rigidly with the phase; the selector shoe state follows the same right-
  continuous branch used by the thermal wall.
- Added focused geometry/energy and reversed-selector tests. All 47 unit tests
  and the scientific validation suite pass.

The named geometry is now available to a renderer, but the Qt preview still
draws particles only. Phase 3 still needs complete-cycle gas/thermal checks.
Phase 4 should next render these components with a fitted camera and add the
apparatus/plot layout.

## Gap closure session 13 — first full apparatus rendering path

- Replaced the particle-only Qt paint path with cylinder, piston, thermal wall,
  piston and selector cams, followers, flywheel/spoke, spring, and load drawing
  from the immutable apparatus snapshot. Triangle orientation marks and the
  selector branch are visible. Generic presets now draw their configured walls.
- Moved the drive schematic below the physical cylinder, leaving the exact
  piston and thermal wall at their collision coordinates. The camera fits a
  full-cycle scene envelope, retains equal world scales on resize, and removes
  the previous hard-coded particle scale.
- Camera tests show all named component paths inside 500×340 and 1050×700
  viewports across several shaft phases for disc and triangle Carnot presets.
  An offscreen PySide6 screenshot from the project virtual environment confirms
  the cylinder, piston, cams, flywheel, spring, and load are drawn at 500×340.
  All 50 unit tests pass in that environment. The full window layout and a
  real display check remain open.

The next Phase 4 slice is the apparatus/plot splitter and worker-owned
simulation; Phase 5 instruments and full-cycle gas/thermal checks remain open.

## Gap closure session 14 — Qt worker ownership and snapshot transport

- Moved physics advancement and stepping out of the GUI timer callback into a
  `SimulationWorker` living in a `QThread`. Play/pause, collision step, and
  fixed-duration advancement are queued signals. The GUI reads only copied
  shape/wall geometry and received snapshots; it no longer accesses the live
  simulation while painting.
- Transport snapshots have bytes-backed read-only NumPy arrays. A frame
  acknowledgement limits queued render frames to one; physics can continue
  while a slower GUI drops intermediate views. Physical advancement remains
  fixed at `0.01` per worker tick, so overload lowers achieved playback speed.
- An offscreen Qt test compares a worker-command Carnot snapshot with a
  headless run at the same seed and physical time. Another test exercises
  pause, resume, explicit step, and clean thread shutdown. All 52 tests pass
  in the project virtual environment.

The next visible Phase 4 slice is the apparatus/plot splitter. Phase 5 still
needs live pressure/area and temperature histories, energy instruments, and
complete-cycle gas/thermal checks.

## Gap closure session 15 — live instrument panel

- Added a roughly two-thirds apparatus / one-third instrument splitter. The
  right panel has a pressure–area trace with points colored by the actual cam
  branch, a translational/rotational temperature history with hot/cold
  reservoir references, key current values, and a scrollable energy and
  numerical-health ledger.
- Pressure is the sum of raw piston-normal impulses divided by piston length
  and a labeled physical-time window (`0.25` maximum). Bounded plot history
  is collected passively in the worker and sent with render frames. Area,
  storage energies, heat, work, and first-law residual come from the same
  simulation state and ledger as collision handling.
- A physical piston-impact test checks the pressure formula and energy split.
  An offscreen Qt test confirms both panels and plots fit at 850×520. The
  reviewed screenshot is `docs/carnot-lab-850x520.png`. All 54 tests pass in
  the project virtual environment.

Remaining: pressure-area work comparison and uncertainty, ideal Carnot
reference, completed-cycle statistics, gas/thermal cycle validation, richer
controls, and the long-run numerical gates.

## Gap closure session 16 — exact completed-cycle ledgers

- The scheduler now records a ledger marker at the initial cam boundary and
  each subsequent completed revolution, before any contact assigned to the
  new right-continuous branch. Markers survive checkpoint/restore and are
  independent of GUI sampling cadence.
- Piston-contact gas energy transfer is accumulated separately from total
  motor/load work. Cycle summaries report signed hot/cold heat, piston work on
  gas, motor work, load output, changes in gas/apparatus storage, and separate
  gas and total first-law residuals. The UI shows completed-cycle count and
  latest-cycle totals; `cycles.json` and `summary.json` export them headlessly.
- Complete controlled cycles with seeded disc (2-particle) and triangle
  (4-particle) gases pass in forward and reversed directions, including hot,
  cold, and piston contacts. Worst cycle residual among those four runs was
  `7.4e-15`, with no penetration above the configured geometry tolerance.
  The exported one-cycle disc run reported one cycle and gas/total residuals
  of `4.44e-16` and `-1.33e-15`. A Qt test also verifies that completed-cycle
  totals reach the instrument panel. All 57 tests pass.

These dilute single cycles are accounting checks, not evidence of quasistatic
Carnot efficiency. Efficiency stays unpublished until transient exclusion,
matched-cycle ensembles, uncertainty handling, and slow-cycle convergence
are validated. Pressure–area work remains a secondary comparison to add.

## Gap closure session 17 — pressure–area cycle work

- For controlled-speed Carnot cycles, integrated the live pressure instrument
  against the exact cam area over each completed cycle. The pressure instrument
  now uses a fixed 0.25 physical-time trailing window, zero-padded at startup,
  so its displayed pressure and work integral use the same definition.
- The work comparison uses raw piston impulses, exact event times, and cycle
  markers. A deterministic bound sums each impulse's smoothing and cycle-end
  clipping effect. The UI labels it as a smoothing bound, not a statistical
  error bar; `cycles.json` includes the comparison values.
- Seed-123 forward/reversed disc and triangle cycles had differences of
  0.2811, 0.1290, 0.1444, and 0.3052 work units, all inside their respective
  bounds of 0.4340, 0.1830, 0.5119, and 0.4690. Narrowing the window from
  0.25 to 0.01 brought the forward disc comparison closer to event work. The
  complete offscreen test suite passes: 58 tests.

Free-shaft comparison awaits stored phase trajectories. The ideal reference,
multi-cycle statistical efficiency, slow-cycle convergence, and remaining
numerical and UI gates remain open.

## Gap closure session 18 — ideal Carnot pressure–area reference

- Added dashed, branch-colored analytical point-gas curves to the live
  pressure–area plot. They use the preset's particle count, hot/cold reservoir
  temperatures, designed cam area profile, and two or three active degrees of
  freedom for smooth discs or rotating triangles.
- The overlay evaluates `P = N T / A` on isotherms and
  `T A^(2/f) = constant` on adiabats. The label shows theoretical
  `ηC = 1 − Tc/Th`; it does not report measured run efficiency. The reference
  is withheld for non-Carnot presets or incompatible cam areas.
- Focused tests verify that all four branches join, the loop integral equals
  the analytical Carnot work, the reference does not change with microscopic
  state or shaft direction, and the Qt pressure–area panel receives it. An
  offscreen 850×520 capture is `docs/carnot-lab-ideal-850x520.png`; the full
  suite passes with 61 tests.

Multi-cycle measured efficiency, uncertainty, slow-cycle convergence, and the
remaining numerical and UI gates remain open.

## Gap closure session 19 — gated completed-cycle efficiency

- Added a signed measured ratio of summed net external output (`load − motor`)
  to summed hot heat over consecutive complete cycles in one direction. The
  default gate skips two transient cycles, requires eight more, and withholds
  a value unless hot heat is safely positive in each of four contiguous cycle
  blocks. Headless and GUI gate parameters are configurable.
- Reports sample count, total storage change, and a four-block delete-one-
  block jackknife standard error. This is descriptive variability, not a
  confidence interval or evidence of quasistatic operation. The UI shows the
  gate state while waiting; `summary.json` replaces the former ungated
  cumulative ratio with the gated result and full report.
- A ten-cycle probe found that accumulated controlled-cam phase timing could
  place the final boundary just after the requested stop. Controlled phase and
  branch event times now anchor to absolute run time. The exact ten-cycle
  regression records ten completed cycles at the requested endpoint with
  first-law residual below `2e-14` in a two-disc seed-123 run.
- Seed-123 ten-cycle runs with two, four, and eight discs correctly withheld
  efficiency because their post-transient hot heat was insufficient or changed
  sign across blocks. Synthetic paired-cycle tests verify a ready estimate,
  transient exclusion, direction resets, heat gating, and export behavior.
  The full offscreen suite passes with 67 tests. A CLI one-cycle export with
  custom gate arguments preserved those settings and wrote null efficiency
  with an `awaiting_matched_cycles` report.

Slow-cycle convergence, larger multi-seed runs, and stationary-state checks
remain open.

## Gap closure session 20 — live diagnostics drawer

- Moved numerical-health readouts into a collapsible instrument drawer and
  added the maximum simultaneous-contact solver residual.
- The worker now measures recent event throughput in events per wall-clock
  second and achieved physical-time playback rate. Rates exclude manual steps
  and display Paused outside continuous playback; the rolling window uses
  recent worker ticks without changing the physical step or event scheduler.
- Offscreen GUI checks cover opening the drawer and the playback-to-pause
  readout transition.

Slow controlled-cycle convergence, multi-seed statistical comparison, and
stationary-state checks remain open.

## Gap closure session 21 — physical shaft speed and reproducible study

- Added positive finite `shaft_speed` to `RunConfig`, `run`, and `gui`. Carnot
  controlled and free shafts use it; `--cycles` now derives nominal physical
  duration from the requested speed instead of a fixed `0.15` rad/s.
- Added `microthermo speed-study` to run controlled Carnot cycles at each
  requested speed and seed. JSON retains exact completed-cycle summaries,
  per-seed efficiency gate, hot heat, external output, storage change,
  numerical health, and the analytical ideal reference. Cross-seed mean and
  standard error appear only with at least two gate-ready seeds.
- A two-disc, ten-cycle probe at speeds `0.075`, `0.15`, and `0.3` rad/s with
  seeds 123–125 completed ten cycles in every run. All nine efficiency reports
  were withheld for insufficient positive hot heat. At the slowest speed,
  post-transient hot heat was negative for all three seeds. The study therefore
  does not establish quasistatic convergence.

Next scientific work needs stationary-state and branch-resolved heat checks
with enough particles and seeds before any comparison to ideal efficiency.

## Gap closure session 22 — branch heat and cycle-energy diagnostics

- The controlled speed study now advances at exact branch boundaries and
  records, per branch and cycle, thermal contact count, signed net heat,
  positive/negative heat contributions, translational/rotational temperatures,
  and gas energy at entry and exit. It reports hot-branch heat reconciliation
  against the existing ledger.
- Added a descriptive post-transient gas-energy screen comparing early and
  late halves of cycle-end energies. A two-standard-error shift is flagged as
  a drift signal; otherwise the result remains `unresolved`. Serially
  correlated cycle samples prevent this screen from certifying stationarity.
- A two-disc probe at `0.15` rad/s reconciled hot heat exactly. Its eight
  post-transient hot branches had 74 thermal contacts, `+61.86` added heat and
  `−84.99` removed heat, yielding `−23.13` net hot heat. Mean temperature on
  entering the hot branch was `3.99` versus hot reservoir `1.5`.
- An eight-disc probe used speeds `0.075` and `0.15` rad/s, seeds 123–124,
  ten cycles, and two transient cycles. All four runs completed ten cycles;
  all four efficiency reports were withheld. Post-transient hot heat ranged
  from `−8.27` to `+10.35`, with positive and negative event heat nearly
  canceling in each run. Mean pre-hot translational temperature was
  `1.61–2.25`, above the `1.5` reservoir. Heat reconciliation was exact to
  printed precision and maximum first-law residual was below `1.3e-13`.
  All four cycle-energy screens were unresolved, not evidence of a stationary
  state. The compact results are in `docs/carnot_branch_probe_session22.csv`;
  reproduce them with `python -m microthermo speed-study --preset carnot_discs
  --speeds 0.075 0.15 --seeds 123 124 --cycles 10 --particles 8`.
- A one-cycle regression confirms exact branch sampling preserves event
  chronology and cycle markers relative to a single advance call.

The hot-heat gate failure is now localized to the thermal exchange data, not
an accounting mismatch. Larger particle ensembles and a stronger stationarity
assessment are still needed before evaluating speed convergence.

## Gap closure session 23 — thermal collision-normal energy

- Added passive thermal event metadata for reservoir temperature and incoming
  and outgoing contact-normal mode energies. The thermal rule samples the
  outgoing mode with expected mean reservoir temperature; the difference of
  outgoing and incoming mode energies equals signed heat, including when a
  boundary moves. The core collision rule and RNG sequence are unchanged.
- The speed study now reports per-branch mode-energy totals and a heat
  reconciliation residual, plus mean incoming/outgoing hot-contact energies
  after transients. Full-cycle forward and reversed regressions check the
  event-level identity.
- Repeated the eight-disc, ten-cycle study at speeds `0.075` and `0.15` with
  seeds 123–124. Across 247–473 post-transient hot contacts per run, mean
  incoming mode energy was `1.456–1.565` and mean outgoing was `1.459–1.571`;
  the hot bath's expected outgoing mean is `1.5`. Their difference accounts
  for each run's small net hot heat, with reconciliation residual below
  `5e-14`. The compact run data are in
  `docs/carnot_normal_mode_probe_session23.csv`.
- Whole-gas translational temperature at hot-branch entry was `1.61–2.25` in
  these runs. It is not interchangeable with the incoming contact-normal
  energy sampled at wall collisions. Near-equal incoming and outgoing means
  explain the hot-heat cancellation; they do not establish thermodynamic
  equilibrium, stationarity, or Carnot convergence.

Next scientific work is to test stationary behavior and heat direction with
larger particle ensembles and longer post-transient windows, then seek
gate-ready multi-seed efficiency estimates.

## Gap closure session 24 — cycle-batch drift screens

- Replaced the earlier per-cycle early/late screen with four contiguous batch
  means in each half of a run. It requires at least 16 cycles after transients
  and compares the early/late shift with a predeclared 10% magnitude tolerance
  and a descriptive two-standard-error margin. It can report `drift_signal`,
  `shift_bounded`, or `inconclusive`; correlated batches mean none of these
  certifies thermodynamic stationarity. The same screen now examines both
  cycle-end gas energy and per-cycle hot heat.
- Ran 18 controlled cycles with 16 discs at `0.15` rad/s for seeds 123 and
  124, excluding the first two cycles. Both runs completed all 18 cycles with
  first-law residuals below `7e-13`. Gas-energy and hot-heat screens were
  inconclusive for both seeds. Their two-standard-error margins exceeded the
  respective 10% shift tolerances. The compact result is
  `docs/carnot_stationarity_probe_session24.csv`.
- Both measured-efficiency reports remained withheld. Seed 123 had `+33.57`
  total post-transient hot heat, but one of its four required contiguous
  hot-heat blocks was `−0.11`. Seed 124 had `−15.76` total hot heat, with
  alternating block signs. This is direct evidence that a positive run total
  alone is inadequate for a stable efficiency estimate.
- The complete offscreen test suite passes with 75 tests.

The present data do not establish a stationary cycle or an approach to ideal
Carnot efficiency. More independent seeds and longer, better mixed runs are
needed; the reference solver's cost now also makes performance work relevant.

## Plan refresh and performance baseline — 2026-09-29

- Replaced `docs/CARNOT_GAP_CLOSURE_PLAN.md` with only outstanding acceptance
  work. Completed sessions remain documented above; the new first priority is
  a performance decision and a physically defined high-particle benchmark.
- Profiled the Python reference backend on an Intel Core i7-11850H. On 32
  discs, the full earliest-event scan and penetration checks dominated. On 16
  triangles, polygon pair TOI dominated. Short unprofiled timings and the
  Numba, algorithmic, and Rust options are in
  `docs/CARNOT_PERFORMANCE_OPTIONS.md`.
- No acceleration backend or particle-size scaling has been selected yet.

## A+B acceleration slice — 2026-09-29

- Selected a fixed occupied area fraction in the Carnot cylinder. Radius is
  `0.025*sqrt(48/N)`; at every count, triangle area occupies about 3.12% and
  disc area about 7.54% of the minimum cylinder area.
- Added a conservative swept circumcircle spatial grid for pair search and
  penetration checks from 32 particles, deterministic pair order, a selectable
  all-pairs reference path, and optional Numba compiled swept-reach and polygon
  transform and witness kernels.
  Cached the immutable cam branch boundaries.
- Installed Numba 0.67.0 and llvmlite 0.49.0 in the project virtual environment,
  pinned the optional Numba range in `pyproject.toml`, exposed `--pair-search`
  and `--numeric-backend`, and added a reproducible short-run benchmark command.
- Grid versus reference tests match event chronology, state, and ledger on 32
  discs and triangles. Randomized swept-grid tests cover 32 and 128 bodies;
  Numba and Python reach masks match. The polygon kernel matches 200 randomized
  triangle pairs in separation, normals, witnesses, and features; its Carnot
  trajectory agrees with the reference within tight tolerances. Full suite:
  81 tests passed.
- Short warm-JIT solver timings and remaining acceleration work are in
  `docs/CARNOT_GAP_CLOSURE_PLAN.md`. The compiled witness kernel makes the
  measured short 32–64 triangle runs about 8× faster than grid Python. Fusing
  transforms and witnesses reduced a short 200-triangle run from 7.06 to 3.72
  wall seconds. A 500-triangle, 0.005-physical-second run took 5.08 wall
  seconds. The 200/500-particle performance target is open.

## A+B acceleration: compiled polygon CCD — 2026-09-29

- Moved the complete polygon pair conservative-advancement loop into an
  optional Numba kernel. It retains the reference touching-contact re-query,
  impact-time bisection, CCD refinement status, and feature data. The Python
  reference remains selectable. Numba uses `fastmath=False`.
- Added direct parity cases for stationary and rotating touch, fast motion,
  and randomized polygon trajectories. Added compiled checkpoint/replay and
  irregular-sampling chronology checks. Backend event order, final state,
  and ledger remain within declared tolerances on the tested cases. The full
  suite passes: 83 tests. A timing-sensitive GUI pause test now waits for the
  worker's paused frame before checking that time stays fixed.
- Warm short probes on the i7-11850H: 200 triangles, 0.02 physical seconds,
  43 events, 2.38 wall seconds (previous compiled witness path 3.72);
  500 triangles, 0.005 physical seconds, 14 events, 2.03 wall seconds
  (previous 5.08). Both reported zero penetration and CCD failures, with
  energy residual magnitude below `1.4e-13`.
- A 200-triangle profile now assigns about 0.10 of 0.21 profiled seconds to
  polygon wall queries and 0.07 to penetration checks on a 0.005-second
  physical interval. Those are the next performance targets. Interactive
  playback and long-run scientific acceptance remain open.

## A+B acceleration: bounded wall search — 2026-09-29

- Added interval-level wall position/speed caches and a conservative swept
  circumcircle rejection before wall TOI. Penetration checks now use circle
  lower bounds before polygon geometry. The `wall_search=all` configuration
  keeps the unpruned reference path.
- Added trajectory parity against unpruned wall search, a high-speed moving
  piston contact, and a deliberately penetrated wall diagnostic. Bounded and
  unpruned paths matched event counts, state, ledger, penetration, and energy
  residuals in the tested short probes.
- On the i7-11850H, 200 triangles/0.02 physical seconds took 2.18 wall
  seconds unpruned versus 1.00 bounded; 500 triangles/0.005 physical seconds
  took 1.99 versus 1.05. The short 200-triangle profile now points to the
  near-wall polygon loop, repeated penetration checks, and grid construction.
  All 85 tests pass. Interactive and long-run scientific gates remain open.

## A+B acceleration: fixed-wall and penetration kernels — 2026-09-29

- Compiled the full polygon CCD loop for walls with constant linear velocity,
  including touching-contact checks and impact-time bisection. Carnot's cam
  piston still uses the Python trajectory path. Added `wall_kernel=python`
  for isolated reference comparisons.
- Batched circumcircle lower bounds for body/wall penetration and compiled
  exact polygon support gaps for near-wall cases. Added
  `penetration_kernel=python` as an independent reference switch.
- Direct fixed-wall tests cover fast motion, rotation, and a linearly moving
  wall. A 48-triangle Carnot trajectory and penetrations near every Carnot
  wall match their Python paths within declared tolerances. All 88 tests pass.
- With other paths fixed, 200 triangles/0.02 physical seconds took 1.02 wall
  seconds with Python wall CCD and 0.81 with the compiled kernel; 500
  triangles/0.005 took 1.10 and 0.99. Compiled penetration checks reduced
  200 triangles from 0.78 to 0.70 seconds and 500 from 0.92 to 0.88.
  These short timings are directional and remain far below interactive rate.
- The next short profile spreads cost across penetration, pair dispatch,
  wall dispatch, and grid construction; there is no single dominant numeric
  kernel left in that sample. Longer parity and scientific acceptance remain
  open.
