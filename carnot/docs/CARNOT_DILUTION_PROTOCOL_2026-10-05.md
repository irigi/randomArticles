# Matched Carnot dilution pilot — predeclared 2026-10-05

The 64-disc pilot doubled count at the preset's fixed occupied fraction and
left all twelve drift screens inconclusive. The preset's 32-disc minimum-area
occupied fraction is about 7.54%; its ideal cam areas use dilute-gas heat
capacity. Finite excluded area is a plausible contributor to the observed
hot-entry excess, but the existing runs do not establish that cause.

`carnot_radius_scale` changes only Carnot particle radius. The default is 1,
so all old configurations retain their geometry. At fixed 32-disc count, use
scales 0.5 and 0.25, corresponding to one quarter and one sixteenth of the
baseline occupied fraction (about 1.89% and 0.47%). Keep both jackets, equal
cam sectors, reservoir temperatures, seeds 123–125, 18 cycles, and speeds
0.075 and 0.15. This is a density perturbation, not an ideal-cam recalibration.
Smaller discs also change collision mixing and initial placement; any response
cannot be assigned uniquely to excluded volume.

Compare each new run with cycles 1–18 of the matching 32-disc both-jacket
control in `runs/gap_followup_2026-10-04/thermo_{slow,fast}_36.json`. Exclude
cycles 1–2. Primary endpoint: paired mean hot-entry translational temperature
change over cycles 3–18. Also report hot heat per disc, thermal contact count,
four-block heat readiness, both fixed drift screens, storage change, collision
health, and seed-level uncertainty. Do not interpret an efficiency ratio as
stationary performance while drift remains inconclusive.

The density hypothesis merits a 36-cycle follow-up at scale 0.25 only if its
mean paired hot-entry reduction is at least 0.10 at **both** speeds, at least
two of three matched seeds improve at **each** speed, all runs finish without
numerical failure, and at least two of three seeds at each speed pass the
positive-hot-heat gate. Scale 0.5 should show an intermediate mean hot-entry
response at both speeds for a monotone density explanation; report a failure
of monotonicity even if scale 0.25 meets the follow-up threshold. If the scale
0.25 threshold fails, stop this variant without a longer unchanged run.

The temporary launcher ran all four independent speed/scale jobs plus the
full unit/scientific validation suite and stored separate stdout, stderr, and
exit codes. It and its batch directory were deleted after review.

## Result and decision

All twelve study runs completed 18 cycles without CCD failure. Maximum
penetration was `7.67e-13`, and the largest absolute final first-law residual
was `1.95e-12`. The four raw JSON records are retained as
`runs/gap_followup_2026-10-04/thermo_dilution_{slow,fast}_{half,quarter}_18_2026-10-05.json`.
All 24 gas-energy and hot-heat drift screens across the four cohorts remained
unbounded: 23 `inconclusive` and one slow-quarter hot-heat `drift_signal`.

| Speed | Radius scale | Ready seeds | Mean hot-entry T | Paired change from baseline ± seed SE | Improved seeds | Mean hot heat per disc |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.075 | 1.0 control | 3/3 | 1.724 | — | — | 4.216 |
| 0.075 | 0.5 | 3/3 | 1.601 | −0.123 ± 0.035 | 3/3 | 5.772 |
| 0.075 | 0.25 | 2/3 | 1.602 | −0.122 ± 0.044 | 3/3 | 5.693 |
| 0.15 | 1.0 control | 3/3 | 1.697 | — | — | 6.241 |
| 0.15 | 0.5 | 3/3 | 1.635 | −0.062 ± 0.050 | 2/3 | 5.859 |
| 0.15 | 0.25 | 2/3 | 1.682 | −0.015 ± 0.029 | 2/3 | 3.415 |

The primary threshold required the quarter-scale mean reduction to reach
0.10 at **both** speeds. It reached 0.122 only at the slow speed and 0.015
at the fast speed. The half-scale response was not intermediate at either
speed: the slow response plateaued, and the fast response reversed. Stop this
dilution variant without a 36-cycle extension.

The saved branches suggest a competing effect. At fast speed, mean compression
temperature rise fell from 0.978 (control) to 0.903 (quarter scale), but mean
cold-branch exit temperature rose from 0.750 to 0.795, leaving hot entry
nearly unchanged. Mean hot contacts per run fell from 8,614 to 7,160. This
supports an inference that dilution changes both compression and thermal
coupling; it does not isolate excluded volume or establish a general density
law. Hot entry remains above the 1.5 bath, and no stationary efficiency trend
is supported.

The independent regression job passed 157 unit tests and all 16 scientific
validation checks. The temporary launcher and batch were removed after
review.
