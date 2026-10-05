# Carnot thermodynamic evidence — 2026-09-30

This study tests heat flow and cycle stability in the controlled-shaft
apparatus. It does **not** establish convergence to ideal Carnot efficiency.

## Method

Each run used 18 complete cycles, seeds shown below, two transient cycles,
and 16 eligible cycles. The existing efficiency gate requires positive net
hot heat in each of four contiguous four-cycle blocks. The gas-energy and
hot-heat drift screens compare four contiguous batch means in each half of
the eligible run. A `shift_bounded` result would still be descriptive, not
proof of stationarity. Shaft speeds are in rad/s. The hot and cold reservoir
temperatures were 1.5 and 0.75, so the analytical reference efficiency is
0.5. No gas temperature clamp was applied. The CLI selected the available
Numba backend (`numeric_backend=auto`) and the default grid pair search.

The complete per-seed cycle, branch, contact, gate, and numerical records are:

- `runs/carnot_thermo_16disc_5seed_3speed_18cycles_2026-09-30.json`
- `runs/carnot_thermo_32disc_3seed_slow_18cycles_2026-09-30.json`
- `runs/carnot_thermo_32triangle_3seed_2speed_18cycles_2026-09-30.json`

The table reports mean total hot heat over the 16 eligible cycles, with
descriptive seed standard error. `Ready` means only that the hot-heat block
gate permits a descriptive efficiency estimate; it does not mean the drift
screens passed.

| Gas | Shaft speed | Seeds | Mean hot heat ± seed SE | Positive totals | Ready | Drift screens |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 16 discs | 0.3 | 5 | −50.48 ± 8.65 | 0/5 | 0/5 | All inconclusive |
| 16 discs | 0.15 | 5 | +0.98 ± 17.36 | 2/5 | 0/5 | All inconclusive |
| 16 discs | 0.075 | 5 | +37.92 ± 16.31 | 4/5 | 0/5 | All inconclusive |
| 32 discs | 0.075 | 3 | +89.45 ± 3.81 | 3/3 | 2/3 | All inconclusive |
| 32 triangles | 0.15 | 3 | −97.79 ± 49.49 | 0/3 | 0/3 | All inconclusive |
| 32 triangles | 0.075 | 3 | −3.06 ± 53.54 | 1/3 | 0/3 | All inconclusive |

## Interpretation

- At 16 discs, each of five paired seeds gained hot heat when the shaft
  slowed from 0.3 to 0.15 and again from 0.15 to 0.075. The paired mean
  increases were 51.46 ± 17.58 and 36.94 ± 12.57, respectively (descriptive
  seed standard errors). This supports a speed effect on hot heat over the
  tested range, while all 15 efficiency estimates remained withheld.
- At 32 triangles, slowing from 0.15 to 0.075 changed mean hot heat by
  +94.73 ± 82.49 across three paired seeds. Two changes were positive and
  one negative. All six triangle efficiency estimates remained withheld.
- At 32 discs and speed 0.075, seeds 124 and 125 passed the hot-heat block
  gate. Their descriptive net efficiencies were 0.183 ± 0.162 and
  −0.026 ± 0.158 by four-block jackknife. The group mean was 0.078 ± 0.105
  by seed standard error, versus ideal 0.5. Seed 123 had positive total hot
  heat but one negative block, so its estimate was withheld. All three gas
  energy and hot-heat drift screens were inconclusive; the two ready seeds
  also changed stored energy by about +10.7 each over eligible cycles.
  These figures cannot support a stationary-engine or Carnot-convergence
  claim.
- Across all 24 runs, the mean translational temperature immediately before
  hot branches exceeded the hot reservoir's 1.5. The thermal wall's mean
  outgoing contact-normal energy stayed near its 1.5 setting; incoming
  contact-normal energy was often similar or greater. This directly explains
  small or negative **net wall heat** in those cases. Overheated gas entering
  the hot sector after compression is a plausible contributing mechanism;
  the study has not isolated cam timing, finite-size effects, and incomplete
  equilibration as separate causes.
- All 24 runs completed 18 cycles, with zero CCD failures. Maximum reported
  penetration was below `9.9e-13`, against a `1e-10` geometry tolerance;
  maximum absolute first-law residual was below `3.2e-12`. Numerical health
  therefore does not explain the failed heat and drift gates in these runs.

## Optional cold-jacket protocol

The baseline cylinder has a cold selector wall at its left edge; its top and
bottom walls are specular. An optional `--cold-jacket` protocol makes those
stationary top and bottom walls thermal at 0.75 only during the cold cam
sector. They are specular in the other sectors. This adds cold contact area
without changing piston motion, the hot boundary, particle count, or the
predeclared heat and drift screens. It is a separate apparatus configuration;
the baseline remains the default. The matched output files are:

- `runs/carnot_cold_jacket_16disc_5seed_2speed_18cycles_2026-09-30.json`
- `runs/carnot_cold_jacket_32disc_3seed_slow_18cycles_2026-09-30.json`
- `runs/carnot_cold_jacket_32triangle_3seed_slow_18cycles_2026-09-30.json`

| Gas | Speed | Baseline → jacket mean eligible hot heat | Heat-gate ready | Drift screens |
| --- | ---: | ---: | ---: | --- |
| 16 discs | 0.15 | +0.98 → +31.94 | 0/5 → 0/5 | All inconclusive |
| 16 discs | 0.075 | +37.92 → +56.75 | 0/5 → 3/5 | All inconclusive |
| 32 discs | 0.075 | +89.45 → +123.28 | 2/3 → 3/3 | All inconclusive |
| 32 triangles | 0.075 | −3.06 → +136.37 | 0/3 → 2/3 | All inconclusive |

At 32 discs and speed 0.075, mean cold-branch exit translational temperatures
for seeds 123–125 fell from 0.829, 0.803, and 0.856 to 0.775, 0.789, and
0.763. Hot heat increased in all three matched seeds. The three jacket runs
passed the four-block positive-hot-heat gate; their descriptive efficiency
mean was 0.337 ± 0.012 by seed standard error, versus the analytical 0.5.
Each run still had inconclusive gas-energy and hot-heat drift screens, and
stored energy changed by −8.4 to −16.8 over the 16 eligible cycles. The
pre-hot gas temperature remained above the hot reservoir in all three runs.
These observations support inadequate cold-sector coupling as one contributor
to the baseline's heat behavior, but do not establish steady operation or
Carnot convergence. The 16-disc jacket cohorts also had inconclusive drift
screens at both speeds; 3/5 slow runs passed the heat gate and none of the
faster runs did.

At 32 triangles and speed 0.075, all three matched seeds gained hot heat:
seed totals changed from −79.97, +99.91, and −29.13 to +88.86, +89.14,
and +231.11. Mean cold-branch exit temperatures fell from 0.851, 0.845,
and 0.847 to 0.782, 0.736, and 0.724. Two jacket seeds passed the heat
gate; their net efficiencies were −0.019 and +0.309. The group descriptive
mean was 0.145 ± 0.164 by seed standard error, too uncertain to support an
efficiency trend. Seed 123 still had a negative hot-heat block. All three
gas-energy and hot-heat drift screens remained inconclusive, and pre-hot
gas temperature remained slightly above the 1.5 hot bath (1.56–1.62).

All 16 jacket runs completed 18 cycles, with zero CCD failures, maximum
reported penetration below `8.1e-13`, and maximum absolute first-law
residual below `5.2e-12`. The matched comparisons make cold-sector coupling
a credible contributor to unstable hot heat; they do not establish stationary
operation, identify a unique cause, or establish Carnot-limit behavior.

## Longer slow triangle jacket observation — 2026-10-01

The matched cold-jacket triangle configuration was extended from 18 to 36
cycles at shaft speed 0.075 for seeds 123–125, keeping two transient cycles,
the same four-block positive-hot-heat gate, and the same drift screens. The
artifact is
`runs/carnot_cold_jacket_32triangle_3seed_speed0p075_36cycles_2026-10-01.json`.
Each seed's first 18 cycle and branch records exactly reproduced the prior
18-cycle artifact.

| Seed | Eligible hot heat (34 cycles) | Heat gate | Gas-energy shift screen | Hot-heat shift screen | Pre-hot translational temperature |
| --- | ---: | --- | --- | --- | ---: |
| 123 | +255.41 | Ready | Inconclusive | Inconclusive | 1.590 |
| 124 | +314.25 | Ready | Shift bounded | Inconclusive | 1.611 |
| 125 | +325.32 | Ready | Inconclusive | Inconclusive | 1.607 |

All four contiguous hot-heat blocks were positive for each seed, improving
heat-gate readiness from 2/3 at 18 cycles to 3/3 at 36 cycles. The block
totals still varied widely: 14.52–132.49 for seed 123, 29.32–137.68 for
seed 124, and 22.85–135.56 for seed 125. All three hot-heat drift screens
remained inconclusive. Only seed 124's gas-energy early–late shift was
bounded under the unchanged descriptive screen; this does not prove
stationarity. Gas temperatures entering hot branches still exceeded the 1.5
hot bath on average. The three descriptive net efficiencies were 0.034,
0.170, and 0.229 (mean 0.144 ± 0.058 seed standard error), versus the ideal
0.5. They are **not stationary-engine efficiency estimates**. Eligible
storage changes were +9.29, −18.24, and −8.39, respectively.

All three runs completed 36 cycles with zero CCD failures and zero reported
penetration; the largest absolute first-law residual was `4.49e-12`.
The study took 2,293 wall seconds on an i7-11850H. Longer observation under
the same protocol improved hot-heat block readiness but did not resolve
hot-heat drift. The next experiment should change one physical feature at a
time, such as thermal contact area or compression timing, while retaining
these screening rules and matched seeds.

### Matched hot-contact-area test — 2026-10-03

An optional `--hot-jacket` makes the top and bottom walls thermal at 1.5
in the hot sector. They retain the cold jacket's 0.75 coupling in the cold
sector and remain specular in both adiabatic sectors. This changes thermal
contact area while holding cam timing, selector wall, particle count, shaft
speed, seeds, and the predeclared heat and drift screens fixed. The new
three-seed, 18-cycle, 16-disc artifact is
`runs/carnot_hot_cold_jacket_16disc_3seed_speed0p075_18cycles_2026-10-03.json`.
Its matched control is the seed-123–125 subset at speed 0.075 of
`runs/carnot_cold_jacket_16disc_5seed_2speed_18cycles_2026-09-30.json`.

| Seed | Cold jacket hot heat | Both jackets hot heat | Paired change | Hot contacts, control → both | Heat gate, control → both |
| --- | ---: | ---: | ---: | ---: | --- |
| 123 | +30.61 | +89.13 | +58.52 | 2,224 → 8,692 | withheld → ready |
| 124 | +91.15 | +98.02 | +6.87 | 2,128 → 8,813 | ready → ready |
| 125 | +94.59 | +47.89 | −46.71 | 2,182 → 8,723 | ready → withheld |

The mean eligible hot heat increased from 72.12 to 78.35, but the paired
changes have mixed signs. Heat-gate readiness stayed at two of three seeds.
All six hot-heat and gas-energy drift screens were inconclusive. In the new
runs, pre-hot translational temperatures were 1.55, 1.95, and 1.95, above
the 1.5 hot bath. Mean incoming hot contact-normal energies were 1.48–1.51,
near the outgoing 1.49–1.52. This is consistent with a large exchange rate
whose positive and negative heat contributions nearly cancel; extra wall
contacts alone did not establish sustained net hot input. The two ready
seeds' descriptive net efficiencies were 0.310 and 0.288, with storage
changes +3.49 and +2.48. No stationary-efficiency claim follows.

All new runs completed 18 cycles with zero CCD failures. Maximum reported
penetration was `7.11e-13`, and maximum absolute first-law residual was
`1.34e-12`. The next controlled protocol change should target compression
timing while retaining this contact-area result as a separate comparison.

### Matched cam-timing test — 2026-10-03

The configurable cam fractions were changed from `(0.25, 0.25, 0.25, 0.25)`
to `(0.25, 0.15, 0.25, 0.35)` in hot, adiabatic expansion, cold, and
adiabatic compression order. This holds hot and cold sector durations, total
cycle time, cam area endpoints, cold jacket, count, speed, seeds, and screens
fixed. It lengthens compression and necessarily shortens expansion; their
effects cannot be isolated from this comparison. The artifact is
`runs/carnot_long_compression_16disc_3seed_speed0p075_18cycles_2026-10-03.json`.
The control is the same slow 16-disc cold-jacket subset used above.

| Seed | Control hot heat | Changed timing hot heat | Paired change | Heat gate, control → changed | Pre-hot temperature, control → changed |
| --- | ---: | ---: | ---: | --- | ---: |
| 123 | +30.61 | +93.07 | +62.47 | withheld → withheld | 1.84 → 1.74 |
| 124 | +91.15 | +63.48 | −27.68 | ready → ready | 1.61 → 1.75 |
| 125 | +94.59 | +37.84 | −56.76 | ready → withheld | 1.75 → 1.97 |

Mean eligible hot heat fell from 72.12 to 64.79, but seed changes have mixed
signs. Heat-gate readiness fell from two seeds to one. All three hot-heat and
all three gas-energy drift screens remain inconclusive. The hot contact
count stayed near 2,100–2,200 per seed; altered heat did not come from a
large contact-area change. Pre-hot temperatures remained above the 1.5 bath.
No stationary-efficiency claim follows.

All three runs completed 18 cycles with zero CCD failures. Maximum reported
penetration was `6.86e-13`, and maximum absolute first-law residual was
`4.55e-13`. Direct compiled and scalar piston collision queries agreed at
phases in all four altered sectors. Independent full-cycle Python/compiled
trajectories diverged in event order for both the default and altered cam,
consistent with amplification of small collision-time differences; this
comparison does not establish long-horizon backend trajectory parity.

### Matched longer cold-sector test — 2026-10-03

The cold sector was lengthened from 25% to 35% of a revolution by shortening
adiabatic expansion from 25% to 15%. Hot exposure, adiabatic compression,
total cycle time, cam area endpoints, cold jacket, particle count, speed,
seeds, and screening rules match the equal-sector control. The effect of
longer cold exposure cannot be separated from the shorter expansion in this
comparison. The 18-cycle artifact is
`runs/carnot_long_cold_16disc_3seed_speed0p075_18cycles_2026-10-03.json`.

| Seed | Control hot heat | Longer cold hot heat | Paired change | Heat gate, control → longer cold | Pre-hot temperature, control → longer cold |
| --- | ---: | ---: | ---: | --- | ---: |
| 123 | +30.61 | +13.84 | −16.76 | withheld → withheld | 1.84 → 1.86 |
| 124 | +91.15 | +71.10 | −20.05 | ready → ready | 1.61 → 1.84 |
| 125 | +94.59 | +58.04 | −36.56 | ready → withheld | 1.75 → 1.83 |

Eligible hot heat fell in all three seeds; its mean fell from 72.12 to
47.66. Heat-gate readiness fell from two seeds to one. All three gas-energy
and all three hot-heat drift screens remained inconclusive. Pre-hot gas
temperature increased in each seed and remained above the 1.5 hot bath.
Mean cold-branch exit temperatures changed from 0.802, 0.746, and 0.754 to
0.793, 0.776, and 0.756; the response was mixed rather than uniformly
cooler. Hot contacts stayed near 2,100–2,300 per seed. These observations
do not support stationary positive hot input or an efficiency comparison.

All three runs completed 18 cycles with zero CCD failures. Maximum reported
penetration was `7.21e-13`, and maximum absolute first-law residual was
`4.22e-13`. A separate one-cycle seed-123 direct run and branch-sampled run
had identical 3,816-event histories and ledgers for this timing.

### Longer slow disc cold-jacket observation — 2026-10-03

The 32-disc cold-jacket cohort at speed 0.075 was extended from 18 to 36
cycles for seeds 123–125, with the same two transient cycles, cam, reservoirs,
and heat and drift screens. The artifact is
`runs/carnot_cold_jacket_32disc_3seed_speed0p075_36cycles_2026-10-03.json`.
Each run's first 18 cycle summaries and 72 branch records exactly reproduce
the earlier matched artifact.

| Seed | Eligible hot heat, 34 cycles | Heat gate | Hot-heat drift | Gas-energy drift | Descriptive efficiency | Storage change |
| --- | ---: | --- | --- | --- | ---: | ---: |
| 123 | +352.66 | Ready | Inconclusive | Inconclusive | 0.353 | −11.83 |
| 124 | +316.28 | Ready | Inconclusive | Inconclusive | 0.354 | −10.61 |
| 125 | +242.80 | Ready | Inconclusive | Inconclusive | 0.270 | −9.56 |

All four contiguous hot-heat blocks are positive in every seed, as at 18
cycles. The descriptive efficiency mean is 0.326 ± 0.028 by seed standard
error, versus the ideal 0.5. These are not stationary-engine efficiency
estimates. Hot-heat early-to-late shifts were +3.88, +3.51, and −1.49 per
cycle; the respective two-standard-error margins were 6.24, 12.85, and
6.46, while the unchanged 10% tolerances were only 1.25, 1.21, and 1.24.
Thus the extra cycles did not bound hot-heat drift. All three gas-energy
screens also remain inconclusive. Mean pre-hot translational temperatures
were 1.75, 1.69, and 1.76, above the hot bath's 1.5.

All three runs completed 36 cycles with zero CCD failures. Maximum reported
penetration was `1.50e-12`, below the `1e-10` geometry tolerance, and maximum
absolute first-law residual was `9.64e-12`. The cohort took about 30 minutes
on the measured machine. Longer observation of this unchanged configuration
strengthens positive-input evidence at one shaft speed but does not establish
stationarity across seeds or a slow-cycle efficiency trend across speeds.

The saved contact totals show why the hot-heat screen remains noisy. Across
the 34 eligible cycles, the three seeds each had about 9,000 hot contacts.
Positive hot-contact heat summed to 6,795–6,875, while negative hot-contact
heat summed to −6,502 to −6,558. The net +243 to +353 is a small difference
between these larger exchanges. Individual eligible-cycle hot heat ranged
from −21.96 to +42.01 across the cohort, with within-seed standard deviations
of 12.3–13.1 per cycle. Mean incoming hot contact-normal energies were
1.459–1.463, close to the 1.5 bath, even while whole-gas pre-hot
translational temperatures averaged 1.69–1.76. This distinction explains
why whole-gas temperature alone does not determine hot-wall heat sign.
It does not establish whether the remaining variability is stationary.

### Matched faster 32-disc cold-jacket cohort — 2026-10-03

To compare shaft speeds at the same particle count, the cold-jacket 32-disc
protocol was measured for 18 cycles at speed 0.15, seeds 123–125. The matched
slow control uses speed 0.075, 18 cycles, the same seeds, apparatus, and
unchanged screens. The new artifact is
`runs/carnot_cold_jacket_32disc_3seed_speed0p15_18cycles_2026-10-03.json`.

| Seed | Fast hot heat | Slow hot heat | Slow − fast | Fast heat gate | Slow heat gate |
| --- | ---: | ---: | ---: | --- | --- |
| 123 | +56.34 | +140.74 | +84.40 | Withheld | Ready |
| 124 | +28.79 | +112.52 | +83.74 | Withheld | Ready |
| 125 | −27.62 | +116.58 | +144.20 | Withheld | Ready |

All three paired hot-heat changes favor the slower speed; the mean change is
+104.11 ± 20.05 by descriptive seed standard error over 16 eligible cycles.
The faster cohort has no heat-gate-ready seed: one seed has negative total hot
heat and the other two have a nonpositive contiguous block. All three faster
gas-energy and hot-heat drift screens are inconclusive, as are those of the
slow control and 36-cycle slow extension. Faster mean pre-hot translational
temperatures were 1.82–1.91, above the hot bath's 1.5. These data strengthen
the speed effect on hot input for this protocol; they do not permit a
stationary-efficiency trend comparison.

All faster runs completed 18 cycles with zero CCD failures, maximum reported
penetration `3.62e-13`, and maximum absolute first-law residual `2.73e-12`.
The run took about 3 minutes 35 seconds. Pressure–area integrals are not
stored in the speed-study artifact; the stationary comparison gate remains
open.

## Remaining thermodynamic work

The existing cold-jacket artifacts also permit a matched particle-count
comparison at speed 0.075, 18 cycles, and seeds 123–125, with identical
screens. At 16 discs, eligible hot heat per particle was 1.91, 5.70, and
5.91; at 32 discs it was 4.40, 3.52, and 3.64. The paired changes have
mixed signs. Heat-gate readiness improved from 2/3 to 3/3, but all six
gas-energy and all six hot-heat drift screens remained inconclusive. Pre-hot
translational temperature remained above the 1.5 bath in all runs. This
comparison was completed from existing artifacts; it does not establish
stationarity or a count-converged efficiency.

The tested changes to hot contact area, compression duration, and cold-sector
duration did not yield reproducible positive hot input with bounded drift. The
36-cycle slow disc extension retains positive hot input but also leaves every
drift screen inconclusive.
The matched 32-disc faster cohort adds a consistent three-seed hot-heat speed
effect, but none of its seeds pass the positive-block gate and its drift
screens are inconclusive.
Then compare net external work and pressure–area work with the analytical
reference as shaft speed falls. Extend the event-work comparison to free-shaft
forward and reversed cycles. The present speed-study artifacts do not include
pressure–area cycle integrals or free-shaft trajectories.
