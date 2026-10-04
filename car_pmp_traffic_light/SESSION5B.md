# Session 5B — Publication figure redesign (capped mechanics)

**Figure-only checkpoint.** Session 5A's independently regenerated capped numerical
results are the exclusive inputs to the *four new named figure pairs*. No paper
paragraphs, tables or citations have been rewritten in this checkpoint.
**`article.tex` and `article.pdf` remain the old, uncapped manuscript.** The old
figure files and generators live under `figures/legacy_uncapped/`, with copies
of the old images at their original paths strictly to preserve compilation of
that unrevised manuscript. Session 6 must replace old `\includegraphics` paths,
figure numbering and captions before the paper can be treated as updated.

## Publication sequence for Session 6

| Sequence | New figure (PNG + vector PDF) | Question addressed |
|---|---|---|
| 1 | `figures/fig1_known28_capped` | How does the known-time driver create the necessary delay? Initial brake pulse, long coast, speed recovery before the 28 s green. |
| 2 | `figures/fig2_information_capped` | How do broad and late priors distribute delay differently under the same realized green at 28 s? Complete red-conditional brown trajectories, contingent green continuations, dashed perfect-knowledge reference. |
| 3 | `figures/fig3_anticipation_capped` | Can a policy accelerate **before any green is possible**, and does the isolated restriction cost anything? Time-domain speed/control and position-domain green branches are in one figure. |
| 4 | `figures/fig4_regret_capped` | When does switch-time knowledge have value? Realized regret over each actual prior support, with the uniform-prior average in each panel. |

## Suggested captions (verify placement when Session 6 builds the paper)

1. **Known-switch control at 28 s.** The capped-mechanics best-found
   deterministic solution brakes almost immediately, coasts for most of the
   approach, then applies throttle beginning around 27.18 s. Free coasting
   would reach the stop line at 25.30 s and therefore cannot alone supply
   the required delay. The trajectory is a feasible numerical local candidate.

2. **Information-dependent allocation of delay.** The broad prior
   `U(0,60)` carries more speed early to preserve the possibility of early
   green but eventually brakes and stops if red persists; the late prior
   `U(26,30)` brakes early and stays rolling. The charcoal dashed trace
   assumes `S=28` is known, solid **brown** is the shared trajectory
   **conditional on continued red**, and **green** branches begin only
   after the green observation. Green at 28 s is emphasized; the late-prior
   plot also has lightly drawn possibilities at 26 and 30 s. The same-S
   comparison is illustrative, not an expected-cost comparison.

3. **Anticipatory acceleration under genuine uncertainty.** For
   `S~U(27,27.5)` with time value 900 CZK/h, the capped model's red-policy
   acceleration begins in the 26.425 s control cell. Requiring no
   *speed-increasing throttle strictly before 27 s* raises expected cost
   from `9.400222554` to `9.410469040` CZK, an advantage of `0.010246486`
   CZK for permitting pre-earliest-green acceleration. The older
   *whole-red-branch restriction* has a different and larger cost difference;
   it is deliberately **not** used for this claim. Dots at 27.03, 27.25,
   and 27.47 s identify realized green observations and their post-green
   colored continuations. These figures demonstrate a resolved local-candidate
   difference, not global optimality.

4. **When switch-time information is worth money.** Conditional realized
   regret is `R(s)=J_feedback(s)-J_known^*(s)`. Its uniform average over
   each prior's actual support gives the capped-model VOIs of
   **0.184982900 CZK** for `U(0,60)` and **0.176970704 CZK** for
   `U(26,30)`. Regret is computed using the multistart best-found 480-cell
   reference and PCHIP interpolation, not a certified globally optimal
   known-time cost curve. It is possible to have positive regret for
   early switches even though the broad policy considers early green.

## Color and plotting contracts

- `#A85C30` (solid brown) ONLY: incomplete-information still-red state path.
- `#147C62` (green): chosen realized post-green continuation. Distinct blue,
  green and purple continuation colors are permissible when several switch
  times appear simultaneously.
- `#344052` (charcoal dashed) in Fig. 2: known switch at 28 s.
- Gray dashed in Fig. 3: counterfactual **no speed-increasing throttle
  before 27 s**, not a no-throttle constraint and not a constraint on the
  entire red period.
- Fig. 4 plots scalar regret, *not a trajectory*, and therefore uses
  neutral charcoal / purple instead of reusing trajectory semantics.
- `model_id=capped_b4_balance`, cap 4.136672736 m/s², brake parameter
  4 m/s²; all green continuations are recomputed using that **same** model.
- Rest is the smooth force-balance numerical convention, not a literal
  claim that a stopped driver applies throttle.

## Source files and reproduction

```
OPENBLAS_NUM_THREADS=1 python make_figures.py
python tests/test_session5b_publication_figures.py
```

`figures/session5b_manifest.json` contains figure source values and numeric
provenance. The figures use these **already verified** datasets:

- `diagnostics/capped_known_sweep/capped_b4_028.00_N480.npz`
- `data/u_0_60_capped.npz`, `data/u_26_30_capped.npz`
- `data/s4_capped900_unrestricted_N1100.npz`
- `data/s4_capped900_before27_N1100.npz`
- `data/voi_capped_regret.npz`, `diagnostics/session5a_voi.json`

The original-model datasets, solver code, VOI scripts, and manuscript remain
unchanged. The next checkpoint is **Session 6: revise LaTeX mathematical
boundary descriptions, model setup, narrative, tables/captions/figure paths,
build the manuscript, and visually verify its PDF**. A systematic fixed-mean
prior-width experiment has not been run; do not claim a general width effect.
