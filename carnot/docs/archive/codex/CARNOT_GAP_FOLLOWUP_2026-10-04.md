# GAP follow-up measurements — 2026-10-04

These runs used the existing code and an i7-11850H with Python 3.12.3. Compact
artifacts and the two replay archives are preserved in
`runs/gap_followup_2026-10-04/`. The temporary launcher, logs, and batch
directory were removed after review.

## Numerical acceptance

Two new 96-triangle, two-cycle default-speed seeds each completed with no CCD
failure or reported penetration. Seed 126 recorded 44,964 events and a final
first-law residual of `-5.40e-13`; seed 127 recorded 46,792 events and
`-9.38e-13`. All eight same-state Python/Numba next-event comparisons matched
for each seed. Direct, sampled, and checkpoint-restored histories matched
exactly at sampling cadences 0.37 and 1.11, respectively. This extends the
multi-seed/cadence evidence but does not cover the remaining hard CCD and
multi-feature cases in the numerical gate.

## Matched hot-and-cold-jacket 32-disc study

Both cohorts used seeds 123–125, 18 cycles, two transient cycles, the existing
positive-hot-heat and drift screens, and the equal-sector cam. Their matched
cold-jacket controls are the existing 32-disc artifacts at speeds 0.075 and
0.15. Adding the hot jacket increased hot contacts about fourfold.

| Speed | Seed | Cold-jacket hot heat | Both-jacket hot heat | Change | Both-jacket heat gate |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0.075 | 123 | 140.74 | 160.01 | +19.27 | Ready |
| 0.075 | 124 | 112.52 | 75.57 | -36.95 | Ready |
| 0.075 | 125 | 116.58 | 169.12 | +52.54 | Ready |
| 0.15 | 123 | 56.34 | 200.21 | +143.87 | Ready |
| 0.15 | 124 | 28.79 | 303.92 | +275.13 | Ready |
| 0.15 | 125 | -27.62 | 95.01 | +122.63 | Ready |

The faster cohort improved from 0/3 to 3/3 heat-gate-ready seeds relative to
its matched cold-jacket control. The slow cohort remained 3/3, but its paired
hot-heat changes had mixed signs. Under the both-jacket protocol, faster heat
exceeded slower heat in two of three paired seeds, reversing the direction of
the earlier cold-jacket speed comparison. **All twelve drift screens** (hot
heat and gas energy for six runs) were inconclusive. Descriptive mean
efficiencies were 0.281 ± 0.093 at speed 0.075 and 0.269 ± 0.073 at 0.15
(seed standard errors), versus the analytical 0.5. These cannot be used as
stationary-engine efficiency estimates or as evidence of a slow-speed trend.
All six runs had zero CCD failures, penetration below `7.7e-13`, and absolute
first-law residual below `4.1e-12`.

## Replay and free shaft

Ten-physical-second, 60-frame/s replay archives completed for 200 and 500
discs. Both contain 601 frames in three chunks. The 200-disc archive has
18,196 events, occupies 13,696,919 bytes, and took about 61 seconds to
write. The 500-disc archive has 67,922 events, occupies 31,588,503 bytes,
and took about 520 seconds. Peak process RSS was about 214 MB and 243 MB,
respectively. These measure storage and offline generation, not sustained
real-display frame rate, seek latency, or interpolation error. The solver
still retains event history while precalculating.

Four 32-disc cold-jacket free-shaft runs, forward and reversed at seeds 123
and 124 for 120 physical seconds, completed with zero CCD failures and
first-law residuals below `1.8e-13`. None completed a full cycle. They
therefore provide no cycle-level pressure-area or efficiency evidence.

## Release test attempt and next action

The initial test command exited before collecting tests because the repository
virtual environment has no `pytest` module. The follow-up used Python's
installed `unittest` runner and the scientific validation CLI, and extended
the both-jacket cohorts to 36 cycles at both speeds. Its results and decision
are in `docs/CARNOT_GAP_FOLLOWUP_36C_2026-10-04.md`. Real-display replay
acceptance and the missing desktop controls remain open.
