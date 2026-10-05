# Next matched physical protocol — predeclared 2026-10-04

This plan uses only the saved 36-cycle, 32-disc, both-jacket records in
`runs/gap_followup_2026-10-04/thermo_{slow,fast}_36.json`. No new simulation
result was examined when setting the comparison and decision rule.

## Diagnosis from saved branches

Across the 34 eligible cycles per run, mean cold-branch exit translational
temperatures are 0.703–0.785, near the 0.75 cold reservoir. The following
adiabatic compression raises them by 0.924–0.995, and mean hot-branch entry
temperatures are 1.636–1.766, above the 1.5 hot reservoir. Hot entry is above
the bath in 23–28 of 34 eligible cycles per run. Mean hot branch exit is
1.462–1.551. Thus the persistent entry excess arises across compression in
these records; its physical causes may include finite-rate work, incomplete
equilibration, finite-size effects, and the chosen area path. Hot-contact net
heat is only 0.5–1.8% of gross absolute exchange, so small changes in the
incoming gas state can change net heat substantially. These observations do
not isolate a cause or imply stationarity.

For the matched 18-cycle comparison, first exclude cycles 1–2 and use cycles
3–18 of the existing 36-cycle controls. Their mean hot-entry temperatures
are, at speed 0.075, 1.684, 1.793, and 1.696 for seeds 123–125, and at speed
0.15, 1.676, 1.621, and 1.794. The control's first 18 cycles exactly match
the earlier 18-cycle artifacts.

## Physical change and comparison

Test both jackets with cam fractions `(0.25, 0.15, 0.25, 0.35)` instead of
`(0.25, 0.25, 0.25, 0.25)`. This gives compression more time at the same
shaft speed while preserving hot and cold exposure durations, all four area
endpoints, count, seed, and reservoir temperatures. Expansion is necessarily
shorter, so any response is attributable to the combined timing change.
The same timing change was previously tested only for 16 discs with the cold
jacket, where hot-heat and hot-entry effects had mixed signs. Its behavior
with 32 discs and both jackets is an open physical comparison, not an assumed
improvement.

Run 18 cycles at each of speeds 0.075 and 0.15, seeds 123–125, 32 discs,
two transient cycles, `efficiency_min_cycles=8`, and both jackets. Compare
each run with the first 18 cycles of its exact matched control. The primary
endpoint is the paired change in mean hot-entry translational temperature
over cycles 3–18. Also report cold-exit temperature, compression temperature
rise, hot branch heat and contacts, positive-hot-heat block readiness, gas
energy and hot-heat drift screens, energy residual, penetration, and CCD
failures. Keep the existing heat and drift screen definitions fixed.

The pilot favors a 36-cycle follow-up only if mean hot-entry temperature
falls by at least 0.10 at **both** speeds and at least two of three paired
seeds improve at **each** speed, without a numerical failure. Otherwise stop
this timing variant and report the negative or mixed result. Even a favorable
pilot cannot establish stationarity: a later 36-cycle comparison must retain
the four-block positive-hot-heat gate and the same gas-energy and hot-heat
drift screens before any stationary efficiency claim. Descriptive seed
uncertainty and storage change remain mandatory.

## Pilot result and decision

The six new runs completed without CCD failures. Raw records are in
`runs/gap_followup_2026-10-04/thermo_compression_{slow,fast}.json`.
The paired endpoint is the cycles 3–18 mean hot-entry temperature difference,
changed protocol minus equal-sector control:

| Speed | Seed 123 | Seed 124 | Seed 125 | Mean paired change |
| ---: | ---: | ---: | ---: | ---: |
| 0.075 | +0.068 | +0.013 | +0.107 | +0.063 |
| 0.15 | −0.017 | +0.048 | −0.064 | −0.011 |

At slow speed, all three hot-entry temperatures rose; at fast speed two fell,
but the mean reduction was only 0.011. Neither speed meets the predeclared
0.10 mean-reduction threshold. Stop this timing variant without a 36-cycle
extension. Cold-exit temperatures and net hot heat also changed in mixed
directions. The slow heat gate fell from 3/3 to 2/3 ready seeds; the fast
gate remained 3/3. All twelve new gas-energy and hot-heat drift screens
were inconclusive. These descriptive heat and efficiency values therefore
do not support a stationary performance claim. The pilot narrows this
specific timing hypothesis but does not isolate the physical source of the
hot-entry excess.
