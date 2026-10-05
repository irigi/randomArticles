# Carnot stationarity diagnosis from saved runs — 2026-10-05

The controlled-shaft, both-jacket 32-disc baseline completed 36 cycles at
shaft speeds 0.075 and 0.15 for seeds 123–125. Each run had positive hot input
in all four uncertainty blocks, but all twelve gas-energy and hot-heat drift
screens were inconclusive. The later 64-disc and dilution pilots were matched
to its first 18 cycles, with cycles 1–2 excluded. Their raw run records and
fixed decision rules are in `docs/CARNOT_PARTICLE_COUNT_PILOT_2026-10-05.md`
and `docs/CARNOT_DILUTION_PROTOCOL_2026-10-05.md`.

| Shaft speed | 32-disc radius scale | Cold exit T | Compression temperature rise | Hot entry T | Hot heat / gross hot exchange | Mean hot contacts per run |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.075 | 1.0 | 0.769 | 0.967 | 1.724 | 0.53% | 17,047 |
| 0.075 | 0.5 | 0.748 | 0.841 | 1.601 | 0.83% | 14,954 |
| 0.075 | 0.25 | 0.781 | 0.847 | 1.602 | 0.85% | 14,341 |
| 0.15 | 1.0 | 0.750 | 0.978 | 1.697 | 1.55% | 8,614 |
| 0.15 | 0.5 | 0.765 | 0.878 | 1.635 | 1.69% | 7,423 |
| 0.15 | 0.25 | 0.795 | 0.903 | 1.682 | 1.02% | 7,160 |

These are three-seed descriptive means over eligible cycles 3–18. The hot
reservoir is at 1.5 and the cold reservoir at 0.75. The ideal cam halves
area during adiabatic compression of spinless discs, for which ideal dilute
gas behavior would double temperature: cold exit at 0.75 would lead to hot
entry at 1.5. The measured baseline enters hotter. Dilution reduces the
compression rise, consistent with a finite-size contribution, but at fast
speed the cold exit becomes warmer as radius decreases. Fewer hot contacts
also accompany dilution. These coupled changes do not isolate a single
physical cause; a systematic numerical error is not proved or excluded by
the small energy residuals alone.

Net hot input is only about 0.5–1.7% of the gross signed hot contact exchange
in these windows. Small changes in incoming state can therefore move net heat
and the efficiency denominator substantially. The 64-disc cohort kept
occupied fraction fixed and did not bound any of its twelve drift screens.
The 32-disc dilution pilot bounded none of its 24 screens and gave one
hot-heat drift signal. Together with the 36-cycle baseline, this supports a
reproducible **failure to establish stationary behavior under the tested
protocols**, not a proof that no stationary Carnot regime exists.

## Next scientific decision

Stop unchanged long runs and further particle-count or radius sweeps. The
next protocol needs to separate cold-sector equilibration from compression
geometry and time: for example, a real cold-sector dwell at fixed area would
increase thermal residence without changing the adiabatic area ratio. Such a
dwell must be represented as actual mechanism motion and logged time, not a
temperature clamp. Before implementing it, define a matched control and a
decision threshold for cold-exit and hot-entry temperatures, positive hot
heat, and the unchanged drift screens. If that physical intervention fails,
retain the experimental label and publish the negative result instead of a
stationary efficiency claim. Pressure–area and free-shaft cycle comparisons
remain downstream of a viable stationary protocol.
