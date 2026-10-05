# GAP follow-up: 36-cycle both-jacket study — 2026-10-04

The saved 36-cycle, 32-disc controlled-shaft studies used hot and cold jackets,
equal cam sectors, seeds 123–125, and shaft speeds 0.075 and 0.15. The first
18 cycles of every run, including branch diagnostics, exactly reproduce the
earlier matched artifacts. The two JSON files containing the six runs are in
`runs/gap_followup_2026-10-04/thermo_{slow,fast}_36.json`.

All six runs completed 36 cycles with zero CCD failures. Maximum reported
penetration was `1.60e-12`; the largest absolute final first-law residual was
`1.29e-11`. Each run passed the existing positive-hot-heat block gate, but
**all twelve drift screens** (gas energy and hot heat for six runs) remain
inconclusive. Positive integrated hot input is therefore established for these
finite windows, not stationary operation.

| Shaft speed | Seed | Post-transient hot heat | Descriptive output / hot heat | Gas drift | Hot-heat drift |
| ---: | ---: | ---: | ---: | --- | --- |
| 0.075 | 123 | 398.55 | 0.389 | Inconclusive | Inconclusive |
| 0.075 | 124 | 273.93 | 0.289 | Inconclusive | Inconclusive |
| 0.075 | 125 | 339.82 | 0.402 | Inconclusive | Inconclusive |
| 0.15 | 123 | 304.19 | 0.244 | Inconclusive | Inconclusive |
| 0.15 | 124 | 499.50 | 0.374 | Inconclusive | Inconclusive |
| 0.15 | 125 | 295.40 | 0.268 | Inconclusive | Inconclusive |

The descriptive seed means are `0.360 ± 0.036` at speed 0.075 and
`0.296 ± 0.040` at 0.15 (seed standard errors), against the analytical Carnot
reference of 0.5. The slower speed has higher hot heat and the higher ratio in
only two of the three matched seeds. The run-to-run uncertainty, nonzero
storage changes, and inconclusive drift screens preclude a stationary
efficiency or quasistatic-limit claim. The late 18 cycles do not consistently
settle: for example, fast seed 123's hot heat falls from 200.21 in the first
16 eligible cycles to 103.98 in the remaining 18, while slow seed 124 rises
from 75.57 to 198.36.

The branch records narrow the physical diagnosis. The preset's hot reservoir
temperature is 1.5, but mean gas temperature at hot-branch entry is
1.64–1.77 across the six runs. Hot contacts add and remove nearly equal gross
heat: net hot heat is only 0.5–1.8% of the sum of positive and negative hot
contact heat. Compression before the hot branch and this large cancellation
make the net input sensitive to fluctuations; these records alone do not
separate cam timing from incomplete thermal equilibration. A useful next
protocol should directly measure and reduce the hot-entry temperature excess
without adjusting the drift gate after the fact.

The release checks in this batch passed: `.venv/bin/python -m unittest discover -s tests -q`
ran 144 tests successfully, and every scientific validation check passed.
The virtual environment lacks pytest, so the built-in unittest runner is the
available full-suite command. These checks do not close the GUI or scientific
acceptance gates.

## Decision and next task

Stop extending this unchanged both-jacket cohort. Use the observed hot-entry
temperature excess to define a matched physical protocol change that separates
compression timing from reservoir coupling, together with a predeclared
stationarity analysis, before another long simulation. The independent replay
delivery task also needs a sustained **real-display** frame-rate, seek, load,
and memory benchmark of the already saved 200/500-disc archives; offscreen
smoke results do not close that gate.
