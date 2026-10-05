# Where the Carnot work stands (2026-10-05)

A condensed digest of the Codex sessions. The original reports are kept, unchanged,
in `docs/archive/codex/`. Design decisions remain in `docs/adr/`.

## What is solid

- **Collision engine.** Event-driven discs and rotating triangles, CCD with
  witness-point contacts, contact clusters, strict failure handling (a severe
  overlap pauses the run and writes a diagnostic checkpoint). Multi-seed runs
  show zero penetration and energy residuals around 1e-12.
- **Bookkeeping.** Heat in from the hot and cold walls, piston work, motor
  work and load output are recorded event by event. Per-cycle first-law
  closure is at floating-point level.
- **Mechanism.** Smooth piston cam with four sectors (hot isotherm, adiabatic
  expansion, cold isotherm, adiabatic compression), thermal selector,
  controlled (motor) and free (flywheel) shaft, reversed cycle.
- **Instruments.** Live pressure–area plot with the analytical ideal-gas
  reference, temperatures, energy plot, cycle ledger, guarded efficiency estimate.
- **Precalculate and replay.** `microthermo precalculate` writes a chunked,
  lossless archive; `microthermo replay` plays it back smoothly with seeking.
  Since archive format v2 (2026-10-05), the replay window shows the same
  instruments as the live window.
- 160 unit tests and 16 scientific checks pass on Ubuntu 24.04 (2026-10-05).

## What is not established

**The simulated engine has not been shown to reach a steady (stationary) cycle.**
Runs of up to 36 cycles with 32 discs, 64-disc runs and diluted-gas runs all
failed the drift checks. Two measured symptoms:

1. After adiabatic compression the gas enters the hot isotherm *hotter* than
   the hot reservoir (≈1.70 instead of 1.5). It first gives heat back to
   the hot wall.
2. Net hot heat is only 0.5–1.7 % of the gross heat that flows back and forth
   at the hot wall, so the efficiency denominator is very noisy.

The 500-triangle, 2-cycle run (`runs/triangle_500_two_cycles_seed123_2026-10-01.json`)
shows the same problem from a different angle:

| Cycle | Q_H (into gas) | Q_C (into gas) | Piston work **on** gas |
| ---: | ---: | ---: | ---: |
| 1 | 347.4 | +16.9 | +109.8 |
| 2 | 49.7 | −0.8 | +72.0 |

In both cycles the gas *absorbs* net work and almost no heat leaves on the cold
isotherm. The gas is still heating up from its initial temperature (1.0,
between the reservoirs at 1.5 and 0.75). It is not yet running as an engine.
A 3-cycle recording that uses the same settings will show this warm-up
transient, not a Carnot engine.

**Update 2026-10-05:** a 96-triangle pilot links the problem to shaft speed.
With jackets, the gas started at T_H and the shaft at 0.0375 rad/s (4× slower
than the default), three consecutive cycles delivered W ≈ 23–26, against an
ideal 21.6, with η ≈ 0.47–0.52, against η_C = 0.5. See `PLAN.md` §3.

## Known rough edges

- `carnot_triangles` still carries an "experimental" label.
- Live 200+ triangle runs are well below real time (500 triangles: about
  0.04–0.05 physical seconds per wall-clock second), so they must be precalculated.
- `run/` in the repository root is a stale output folder from 2026-09-29.
- The README is long and was written as a release-evidence document.
