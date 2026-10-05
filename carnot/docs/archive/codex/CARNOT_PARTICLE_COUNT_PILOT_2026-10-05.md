# Particle-count thermodynamic pilot — predeclared 2026-10-05

The 32-disc, both-jacket, equal-sector cohorts at shaft speeds 0.075 and 0.15
have positive hot heat in all four blocks for all three seeds through 36 cycles,
but all twelve gas-energy and hot-heat drift screens are inconclusive. Hot-entry
temperature exceeds the 1.5 reservoir mean, and hot contact heat nearly
cancels between gains and losses. Changing compression time did not reduce the
hot-entry excess at both speeds. Repeating those protocols longer has no clear
decision value.

## Matched change

Run 64 discs with both jackets and equal cam sectors for 18 cycles at each of
0.075 and 0.15 shaft speed, using seeds 123–125. Keep two transient cycles,
the existing four-block positive-hot-heat gate, and the existing 10% gas-energy
and hot-heat drift screens. The comparison is to cycles 1–18 of each matched
32-disc cohort in `runs/gap_followup_2026-10-04/thermo_{slow,fast}_36.json`.
Particle radius scales with count in the preset, so occupied area fraction
remains fixed. The pilot asks whether finite-count fluctuations dominate the
inconclusive drift screens and hot-entry excess; a difference is not assumed.

## Analysis and decision rule

For each speed and seed, report completion, numerical health, efficiency-gate
status, both drift statuses and their mean shifts and margins, mean hot-entry
translational temperature over cycles 3–18, hot heat per particle, and
post-transient storage change per particle. Compare paired 64- versus 32-disc
changes and report seed means with standard errors. Preserve the same screen
definitions; do not infer stationarity from an efficiency-gate `ready` result.

Consider a 36-cycle 64-disc follow-up only if all six runs finish without a
numerical failure, at least two of three seeds at **each** speed pass the
positive-hot-heat block gate, and at least four of the twelve drift screens
say `shift_bounded` at 18 cycles. Otherwise stop the particle-count variant
and document the tested model's failure to establish a stationary efficiency
trend. A favorable pilot still requires longer-run drift checks, pressure-area
work, and free-shaft comparisons before any stationary performance claim.

The temporary launcher and its output directory were deleted after the
results were condensed and reviewed. The batch also ran an independent
fresh-venv install of the optional GUI
and accelerator dependencies on the Ubuntu 24.04 host, then runs the full
unittest and scientific suites. Its package list, stdout, stderr, and exit
code provide release evidence; a fresh venv is narrower than a clean OS image.

## Result and decision

All six 64-disc runs completed 18 cycles and passed the four-block positive
hot-heat gate. There were zero CCD failures, maximum penetration was
`8.26e-13`, and the largest absolute final first-law residual was `4.27e-12`.
The retained raw results are
`runs/gap_followup_2026-10-04/thermo_64disc_{slow,fast}_18_2026-10-05.json`.

| Speed | Ready seeds | Bounded drift screens | Mean hot-entry T, 64 discs | Paired change from 32 discs | Hot heat per disc, 64 / 32 | Descriptive efficiency, 64 discs |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.075 | 3/3 | 0/6 | 1.672 | −0.052 ± 0.023 | 6.18 / 4.22 | 0.409 ± 0.021 |
| 0.15 | 3/3 | 0/6 | 1.722 | +0.025 ± 0.094 | 4.74 / 6.24 | 0.358 ± 0.020 |

The paired-change uncertainty and efficiency uncertainty are seed standard
errors for only three seeds. Hot-entry temperatures remain above the hot bath
temperature of 1.5 at both speeds. Hot heat per disc increases at slow speed
but decreases at fast speed, and no gas-energy or hot-heat drift screen is
`shift_bounded`. Every drift screen is `inconclusive`; none signals established
drift either. The paired post-transient storage change per disc is
`−0.242 ± 0.280` at slow speed and `−0.770 ± 0.309` at fast speed, again
descriptive. The six runs contain 1,940,890 collision events in total.

The predeclared threshold required at least four of twelve bounded drift
screens; the observed count was zero. **Stop this particle-count variant**
without a 36-cycle extension. The 64-disc cohort reinforces positive hot
input over the measured window but does not establish stationary operation or
an efficiency trend toward the 0.5 analytical reference. Further unchanged
long runs are not justified by this pilot.

The independent fresh-venv job passed on Ubuntu 24.04.5 LTS, Python 3.12.3:
NumPy 2.4.6, Numba 0.67.0, llvmlite 0.49.0, PySide6 6.11.2, and PyQtGraph
0.14.0 installed; the headless import did not load Qt; 156 unit tests passed;
and all 16 scientific validation checks passed. The GUI tests ran with the
offscreen Qt platform. This is a fresh environment on an existing host, not
a clean operating-system image or a second real-display acceptance run.
