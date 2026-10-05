# Plan — Carnot visual redesign and 500-triangle replay

Decisions (2026-10-05): physics option **C** (gas starts at T_H, hot and cold
jackets), cycle dial **yes**, temperature-vs-time strip **kept**.
The temperature ratio stays at 2 (see §3).

## 1. New layout — done

```
┌──────────────── left ────────────────────────────────┬── right ─────────────┐
│  cylinder + gas + piston (full barrel dashed)  │ dial │ STATE                │
├──────────────────────────┬───────────────────────────┤ CUMULATIVE           │
│   p–V (measured, by      │   T–S (ideal-gas entropy  │ LAST CYCLE           │
│   branch; ideal dashed)  │   estimate; ideal rect.)  │ EFFICIENCY           │
├──────────────────────────┴───────────────────────────┤ T(t) strip           │
│  current cycle Q_H Q_C W  │ table: one row per cycle │ ▸ Mechanism          │
│                           │                          │ ▸ Numerical diag.    │
└──────────────────────────────────────────────────────┴──────────────────────┘
```

- The flywheel, cams, selector, load and spring drawings are removed. The cylinder is drawn
  alone; active thermal jackets colour its top and bottom walls.
- The cycle dial shows the four cam sectors, the shaft phase, and progress
  through the current sector.
- p–V and T–S show the previous and current cycle in branch colours and older
  cycles in dim grey. Pressure samples are hidden until the 0.25 s window
  first fills.
- T–S uses `S/N = ln A + (f/2) ln T`, with `T = 2E/(fN)` taken from the gas energy. This is
  an estimate, not a measurement; S = 0 is the start of the hot isotherm.
- Scoreboard sign convention: W > 0 means work delivered by the gas. Warm-up
  cycles (the efficiency gate's transient cycles) are greyed.
- Right column: 17 rows in four groups, plus collapsed Mechanism and
  Diagnostics sections. The inspector appears only when something is selected.
  At 1280×800 it fits without scrolling (tested).
- Plots redraw at most about 7 times a second while playing. The particles,
  dial and numbers update every frame. Paused or stepped frames always redraw.
- Non-Carnot presets keep the cylinder-free scene and hide the diagrams.

Code: `ui/lab_view.py` (`LabView`, `CycleDial`), `ui/instruments.py`
(`CycleDiagrams`, `ReadoutPanel`, `InstrumentSeries`).

## 2. Replay shows everything — done

- Archive **format v2**: per-frame instrument columns (`i_*`) plus the cycle
  markers and pressure–area comparisons in `manifest.json → instruments`.
  Version-1 archives still open, showing only the cylinder and a one-line readout.
- `LiveInstruments.observe(..., events=...)` accepts drained events, so
  precalculation records exactly what the live window shows. A test checks
  this field by field, including the cycle ∮P dA integral.
- The replay window uses the same `LabView`. Long replays are thinned to
  3000 plot points; the most recent 1500 stay at full resolution.
- `precalculate` gains `--cycles`, `--initial-temperature`,
  `--temperature-ratio` and `--pair-search` (default `sweep`, 2.1× faster for
  500 triangles), and prints progress to stderr after each chunk.
- `gui` gains `--seed --particles --hot-jacket --cold-jacket
  --initial-temperature --temperature-ratio`.

## 3. Physics: what makes it run as an engine

Pilot: 96 triangles, seed 123, 3–4 cycles each. Raw results are in
`docs/results/carnot_pilot_96tri_2026-10-05.jsonl`; the script is
`examples/carnot_cycle_pilot.py`. The ideal values are W = N·ΔT·ln(A2/A1) = 21.6 per
cycle and η_C = 0.5.

| Variant | Shaft (rad/s) | W by gas per cycle | Verdict |
|---|---:|---|---|
| A: default | 0.15 | −9.4, −14.7, +1.0, −34.8 | absorbs work |
| C: jackets + start at T_H | 0.15 | +17.5, −17.7, +5.1, −8.7 | erratic |
| C | 0.075 | +1.6, +7.0, +0.3, +2.3 | weak engine |
| **C** | **0.0375** | **+22.7, +24.3, +25.7** (η 0.52, 0.51, 0.47) | **Carnot-like** |
| C, T_H/T_C = 3 | 0.15 | −110, −158, −66, −112 | gas overheats (T → 2.5) |
| C, T_H/T_C = 4 | 0.15 | −290, −275, −326, −340 | gas overheats (T → 3.4) |

Interpretation: at 0.15 rad/s the piston moves at roughly 0.3× the thermal
speed, so compression is not quasi-static and heats the gas beyond the
adiabat. The gas then enters the hot isotherm above T_H (the "hot entry"
symptom in `STATUS.md`). A larger temperature ratio needs a larger
adiabatic area ratio (A3/A2 = ratio^(f/2)), which amplifies the same error.
Slowing the shaft fixes it. A 4× slower shaft gives near-ideal work and
efficiency in all three cycles. This is a 96-particle, single-seed pilot,
not a statistical claim.

## 4. Next ideas

1. Record the default setup (A) as a comparison archive. Watching the two
   side by side is the clearest lesson in why a Carnot engine must run slowly.
2. A speed sweep at 0.0375 with more seeds and cycles, to turn §3 into a
   real result, and to drop the "not stationary" caveat if it holds.
3. Make 0.0375 rad/s with jackets the default for `carnot_*` presets.
4. Export a replay to MP4 (offscreen render piped to ffmpeg).
5. Briefly flash particles leaving a hot or cold wall, to show where heat enters.
6. Shorten the README into install, main commands and the replay workflow.
