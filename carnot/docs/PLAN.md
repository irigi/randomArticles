# Plan — Carnot visual redesign and 500-triangle replay

Status: **draft for discussion**. Items marked ❓ need your decision.

## 1. New layout of the Carnot window

### 1.1 What goes away

From the apparatus drawing, remove the flywheel and its spoke, the piston cam,
cam follower, selector cam and shoe, the load arm, the torsion-spring marker,
and their text labels. The cylinder, gas, piston and coloured thermal walls
stay. The mechanism still runs underneath. Its energies stay in the ledger,
and clicking the cylinder or piston still opens the inspector.

❓ Add a small **cycle dial** to the corner of the cylinder view as a replacement
(optional): a ring cut into four coloured sectors (hot, expansion, cold,
compression) with a needle at the current shaft phase. It shows the branch and
how far through it the cycle is, which the cams showed only indirectly.

### 1.2 Proposed arrangement

```
┌─ toolbar ──────────────────────────────────────────────────────────────────┐
├──────────────── left: ~65 % ─────────────────────────┬── right: ~35 % ──────┤
│                                                      │ STATE                │
│        cylinder + gas + piston      (◔ cycle dial)   │ CUMULATIVE LEDGER    │
│                                                      │ LAST CYCLE           │
├──────────────────────────┬───────────────────────────┤ EFFICIENCY           │
│   p–V  (measured +       │   T–S  (gas estimate +    │                      │
│   ideal reference)       │   ideal rectangle)        │ T(t) strip plot      │
├──────────────────────────┴───────────────────────────┤                      │
│  CYCLE SCOREBOARD   W_net  Q_H  Q_C  η  vs  η_Carnot  │ ▸ Mechanism          │
│  (current cycle running + table of completed cycles) │ ▸ Numerical diag.    │
└──────────────────────────────────────────────────────┴──────────────────────┘
```

- **p–V**: the existing measured trace, coloured by branch, plus the dashed
  ideal reference. With much more room it becomes the main plot.
- **T–S** (new): x axis is the gas entropy estimated from the measured
  temperature and area, `S/N = ln A + (f/2) ln T` (k_B = 1, additive
  constant dropped; f = 2 for discs, 3 for triangles). The ideal Carnot
  rectangle is drawn dashed. This is an **ideal-gas estimate** derived
  from measurements. It is not a measured entropy, and the axis label will
  say so. Its strength is that a lagging or overshooting temperature shows up
  directly as a non-rectangular loop, which is the problem in the STATUS document.
- **Cycle scoreboard** (new): large numbers for the cycle in progress
  (Q_H, Q_C, W by gas so far), and a compact table of completed cycles:
  `# | Q_H | Q_C | W_net | η = W/Q_H | η_C`. Engine-sign convention:
  positive W means work *delivered*, so a cycle that consumes work turns red.

### 1.3 The right-hand column, trimmed to one screen

Today's column has a summary line, the inspector, three plots and 27 rows,
plus 6 diagnostics. Proposal:

**Keep, grouped (16 rows)**

| Group | Rows |
|---|---|
| State | Branch · phase · time · cycles completed |
| | Area A · Pressure P (window) |
| | T translation · T rotation |
| | Reservoirs T_H / T_C |
| | Gas energy |
| Cumulative | Q_H · Q_C |
| | Motor work · Load output |
| | First-law residual |
| Last cycle | Q_H / Q_C |
| | Work by gas: event ledger vs ∮P dA (± bound) |
| | External output |
| | Total first-law residual |
| Efficiency | η_net ± block SE  vs  η_Carnot |
| | Gate status (e.g. "3/8 cycles") |

**Remove or move**

| Row | Why it is least relevant | Goes to |
|---|---|---|
| Summary line at top | Duplicates the rows below | removed |
| Piston energy, Shaft inertia energy, Spring energy | Near zero, or set by the motor, in controlled mode | collapsed "Mechanism" section |
| P dA − event work | The difference of two rows that are both shown | merged into the "Work by gas" row |
| Last cycle gas residual | The total residual is the one that matters | Numerical diagnostics |
| Efficiency cycles / storage ΔE | Explains the gate, so it belongs in a tooltip | tooltip on the η row |
| Events | Already in the status bar | status bar only |
| Inspector box | Takes space even when nothing is selected | appears only on selection |
| Energy & ledger plot | Shows the same numbers as the ledger and the scoreboard | removed (still exported) |
| p–V plot | — | moves to the left |

❓ Keep the **temperature-vs-time** strip on the right? I recommend keeping it.
It is the clearest view of "the gas never reaches the reservoir temperature",
and the T–S diagram does not show time.

### 1.4 Code changes (high level)

- `ui/instruments.py` → split into `CycleDiagrams` (p–V, T–S, scoreboard; left)
  and `ReadoutPanel` (right column).
- `ui/main_window.py` → `ApparatusView` gains a `minimal` mode that draws only
  the cylinder; a new `CarnotLabView` widget combines the apparatus view,
  diagrams and readout. The same widget will be used by **both** the live window and the replay window.
- Non-Carnot presets (gas box, labyrinth) keep the current layout.
- Update the `tests/test_ui_render.py` checks; add tests for the entropy estimate and the scoreboard sign convention.

## 2. Replay shows everything the live window shows

**Problem.** The archive stores particle snapshots and four cumulative ledger
numbers. It does not store the pressure window, area, cycle markers, cycle
summaries or efficiency report, so the replay window cannot draw the
instruments.

**Fix.** While precalculating, run the same `LiveInstruments` observer on every
saved frame and store its output:

- Per frame, new columns in each chunk: area, pressure, motor work,
  piston/flywheel/spring energy, completed cycles, and diagnostics.
- Once per run, a new `cycles.json` file with every cycle marker and its
  pressure–area comparison. The replay rebuilds the "last cycle" and efficiency
  rows from it at any seek position.
- Bump the archive format to version 2. Version-1 archives still open, with the
  old one-line readout.

One trap to handle: `LiveInstruments.observe()` reads `simulation.events`, but
precalculation drains that list after every frame. The observer needs an
explicit event feed instead of a cursor into that list.

The replay window then becomes the same `CarnotLabView` with a time slider. The
plots show history up to the current replay time, so seeking backwards works.

**Cost estimate for 500 triangles × 3 cycles** (shaft 0.15 rad/s → 125.7 physical s):
at the measured 0.04–0.05 physical s per wall-clock second, recording takes
**about 45–55 minutes**. The archive will be roughly 250–350 MB at 60 frames/s;
`--fps 30` would halve the frame part of that.

## 3. ❓ Physics choice for the showcase recording

As shown in `STATUS.md`, with the default settings the 500-triangle gas spends
the first cycles heating up and *absorbs* work. Three cycles will look like
"not an engine". Options, from least to most change:

| Option | What changes | Expected effect |
|---|---|---|
| A. As is | nothing | Honest, but shows the warm-up transient |
| B. Start the gas at T_H | new `--initial-temperature` flag (today one `temperature` sets both the gas and the reservoirs at 1.5× / 0.75×) | Removes most of the cycle-1 heating transient |
| C. B + thermal jackets | `--hot-jacket --cold-jacket` make the top and bottom walls thermal during their sector | More wall contact, so the gas follows the reservoir temperature better on the isotherms |
| D. C + slower shaft | `--shaft-speed 0.075` | Doubles the time on each isotherm. Recording takes about 2 hours |

My recommendation: record **C**, and later run **A** for comparison. Both recordings are
useful, and the difference between them is itself instructive.

## 4. Other improvements you might like

1. **Rate-of-heat colouring**: briefly flash each particle orange or blue when it
   leaves a hot or cold wall. You would see where and when heat enters.
2. **Export a replay to MP4** (offscreen render → ffmpeg) for sharing,
   instead of the current PNG-frame screen recording.
3. **Shorter README**: move the long evidence sections into `docs/` and keep
   install, the main commands and the replay workflow.
4. **Cold-dwell experiment** (Codex's last open idea): give the cold isotherm a
   longer sector (`--cam-fractions`) so the gas has time to cool before
   compression. This targets the "hot entry above T_H" symptom directly.

## 5. Order of work

1. Item 2 (archive v2 and the instrument feed) first, so you can start the
   ~1 h recording while the UI is being built.
2. Item 1 (layout), developed against short live runs and the new archive.
3. Hook the replay window up to the shared view, then run the full replay check.
4. Update the README with the final commands.
