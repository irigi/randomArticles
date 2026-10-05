# Carnot release evidence map — 2026-10-05

This maps the Carnot-specific requirements in sections 9, 13 and 14 of
`microscopic_thermodynamics_implementation.md` to reviewable evidence. A test
of accounting or UI behavior does not by itself establish a stationary heat
engine. Other experiments in the brief need separate release mapping.

| Brief requirement | Evidence | Status or limit |
| --- | --- | --- |
| §9.1 wall impulse/pressure, local heating, triangle equipartition introductions | `README.md` presets; `tests/test_live_instruments.py`, scientific validation | Implemented; inspect lesson wording during final release review. |
| §9.1 2D disc/triangle ideal-gas formulas and finite-size caveat | `tests/test_ideal_reference.py`; `docs/CARNOT_THERMODYNAMIC_EVIDENCE_2026-09-30.md` | Analytical reference shown separately from measured data. |
| §9.2 piston/selector cams and four physical branches | `tests/test_apparatus.py`, `tests/test_carnot_cycle.py`, `tests/test_ui_render.py`; real-display acceptance in `docs/CARNOT_GAP_CLOSURE_PLAN.md` §4 | Tested controlled and reversed selection. |
| §9.2 area calibration by gas degrees of freedom, smooth periodic cam, unequal sector durations | `tests/test_ideal_reference.py`, `tests/test_speed_study.py`; `src/microthermo/core/mechanisms.py` | Preset ratio and reservoir temperatures are fixed multiples of the configured initial temperature. |
| §9.2 physical switching without temperature clamp | `tests/test_speed_study.py`; branch/contact diagnostics in the retained thermodynamic JSON | Finite-speed hot entry exceeds reservoir temperature in measured cohorts. |
| §9.3 controlled and free shaft, motor/load/spring bookkeeping | `tests/test_carnot_cycle.py`, `tests/test_apparatus.py`; free-shaft probe in `docs/CARNOT_GAP_FOLLOWUP_2026-10-04.md` | Mechanical contracts tested; sampled free-shaft gas runs did not complete a cycle. |
| §9.4 pressure–area, temperature, component-energy and reservoir instruments | `tests/test_live_instruments.py`, `tests/test_ui_render.py`; `runs/gap_followup_2026-10-04/desktop_files_acceptance.png` | Real-display check covers the measured xcb desktop. |
| §9.4 signed heat, event/mechanism work, first-law closure, pressure–area integral as secondary estimate | `tests/test_carnot_cycle.py`, `tests/test_efficiency.py`; `src/microthermo/measurements/pressure_area.py`; retained pilot JSON | Cycle accounting passes; ensemble pressure–area comparison remains open. |
| §9.4 guarded efficiency, ideal efficiency, uncertainty, storage change | `tests/test_efficiency.py`, `tests/test_speed_study.py`; `docs/CARNOT_PARTICLE_COUNT_PILOT_2026-10-05.md` | Descriptive ratios only; all 64-disc drift screens inconclusive. |
| §9.4 reservoir entropy and nonclosed-run qualification | `src/microthermo/measurements/ledger.py`; `README.md` | Verify user-facing distinction before release. |
| §9.5 shaft-speed, gas shape, particle-count and density comparisons | `docs/CARNOT_THERMODYNAMIC_EVIDENCE_2026-09-30.md`; `docs/CARNOT_PARTICLE_COUNT_PILOT_2026-10-05.md`; `docs/CARNOT_DILUTION_PROTOCOL_2026-10-05.md` | Multi-seed windows measured; stationary slow-cycle trend unestablished. |
| §9.5 reservoir ratio, load/inertia/spring, reversed refrigeration comparisons | `tests/test_carnot_cycle.py`, `tests/test_ideal_reference.py`; `README.md` reversed example | Component behavior tested; the requested systematic multi-seed comparisons and coefficient-of-performance evidence remain open. |
| §9.5 default 200–500 triangles, separate benchmarks | `docs/CARNOT_GAP_CLOSURE_PLAN.md` §§1, 3, 6; retained performance/replay artifacts | Measured numerical and replay gates pass; live 200-triangle target remains deferred and below target. |
| §13.1 Ubuntu install and dependency compatibility | Fresh venv on Ubuntu 24.04.5, Python 3.12.3; details in `docs/CARNOT_PARTICLE_COUNT_PILOT_2026-10-05.md`; clean-image protocol in `docs/CARNOT_CLEAN_UBUNTU_PROTOCOL_2026-10-05.md` | 156 tests and 16 scientific checks passed in the fresh venv; clean OS image still untested. |
| §13.3 throughput and responsiveness | `docs/CARNOT_GAP_CLOSURE_PLAN.md` §§1, 4, 6; retained standard benchmark and desktop screenshots | Measured cases only; aspirational full benchmark matrix remains incomplete. |
| §14.1–14.2 geometry, contact, conservation, checkpoint validation | `docs/CARNOT_GAP_CLOSURE_PLAN.md` §3; `tests/test_carnot_cycle.py`, `tests/test_backend_branch_parity.py` | Measured numerical gate passes; no million-impact isolated trajectory claimed. |
| §14.3 statistical mechanics and Carnot slow-cycle convergence | `docs/CARNOT_THERMODYNAMIC_EVIDENCE_2026-09-30.md`; `docs/CARNOT_PARTICLE_COUNT_PILOT_2026-10-05.md` | Statistical validation exists; stationary Carnot convergence remains unestablished. |

Release acceptance remains open until the untested comparisons and scientific
claim are resolved or explicitly scoped out in the shipped product, and a
clean operating-system installation plus the rest of the brief's requirement
map are reviewed. The `carnot_triangles` experimental label remains accurate.
