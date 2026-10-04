# Checkpoint 3 — known-time reference, VOI and information regret

Third session of six. Builds on Session 2 (`SESSION2.md`), preserving its
verified uncertain-policy arrays and multistart/feasibility-tested solvers.

**This is still a numerical checkpoint, not a revised publication**:
`article.tex`, `article.pdf`, all `figures/` publication images, the
anticipatory counterfactuals and the mechanical-sensitivity arrays are carried
forward for Sessions 4–6. In particular, the PDF still contains the **old**
Table II numbers, which are now superseded numerically in `results.json` and
`diagnostics/session3_voi.json`.

## Main numerical estimates, original uncapped mechanics, 300 CZK/h

| Prior | Best local known-S average | Feedback expected | VOI | old VOI |
|:--|--:|--:|--:|--:|
| U(0,60) | 4.275751587 | 4.407832740 | **0.132081152** | 0.129 |
| U(26,30) | 3.671992587 | 3.843013166 | **0.171020579** | 0.171 |

All costs are CZK. The known-S solution is based on a grid at every integer
switch time from 22–60 s, plus 26.5/27.5/28.5/29.5 s (each at 480 cells),
and the exact free-flow threshold at D0/VMAX = 21.6 s. Known costs below this
threshold equal the free-flow route cost. The curve between known-S knots is
integrated using a *shape-preserving piecewise cubic*, not guessed from the
nearest sampled policy. It remains an approximation to a local-best solution,
not a global-optimality certificate. Distinguish these numeric estimates from
the publication's earlier 336-cell S=28 illustrative trajectory.

The previous `value_of_information.py` had two errors: it divided a cost
integral over 26–60 s by four for the 26–30 s prior, and described 0.5 s
sampling while computing at a 1 s step. The script now integrates each prior
over its **exact support**, keeps explicit known-time knots and rejects missing
known-solution files rather than silently solving single-start problems.

## Conditional regret rather than speculation about stops

At a realized green time s define

    R(s) = J_feedback(s) - J_known_best(s).

The feedback realization is reconstructed from the respective Session 2
red-conditional policy and green value: each transcription-cell running cost
is integrated up to s; the state is interpolated within the switch cell.
Eight-point Gauss-Legendre integration over every red-branch cell gives:

- Broad reconstruction minus saved feedback expected cost: -0.00000609 CZK.
- Late reconstruction minus saved feedback expected cost: -0.0000000323 CZK.

The near agreement independently checks the policy expectation calculation,
and numerical regret stays nonnegative at all saved grid samples within the
interpolation tolerance. Integrated regret and `feedback_NLP - known` differ
only by the small quadrature mismatch above.

Broad-prior regret peaks near S=26 s at about 0.431 CZK and is already positive
for earlier green times (e.g. roughly 0.053 CZK at S=10 and 0.207 CZK at
S=20). In the late-prior window regret decreases from about 0.340 CZK at
26 s toward 0.012 CZK at 30 s. These statements are descriptive of the
**computed** policies and reference curve; they do not assume that late
known-S switches *force* stop/wait. Indeed, S=40 and S=50 now both retain
rolling deterministic candidates.

## Reliability, resolution and model limitations

- Multi-start physical initializations and selected reverse continuation passes
  were used on the 480-cell 22–60 s grid. Every reference point is saved
  individually under `diagnostics/known_sweep/` and is resumable.
- At S=40, N=480 gives 4.95568147 while an N=960 warm-started check gives
  4.95562795 CZK; at S=50, 5.91244709 versus 5.91233961 CZK.
- At S=57, N=480 gives 6.55633279 versus 6.55634699 at N=960.
- At S=60, the N=480 candidate is 6.82381257 CZK while one N=960 candidate
  is 6.82835959 CZK. Near the near-rest restart the *uncapped* low-speed
  acceleration law and time discretization make some solutions sensitive to
  grid size, candidate basin and whether a trajectory momentarily stops. As
  an especially stark example, N=240 produces a lower candidate at S=60
  that is not reproduced when refined to N=480; it must not be promoted to
  a proven lower physical optimum. Fully rigorous convergence/globality is
  **not** claimed. Session 4 revisits mechanical realism, which will require
  a new main-model VOI calculation if the baseline mechanics are changed.
- The high-detail S=28 numerical example from the paper is not silently
  replaced: the 480-cell VOI reference and the 336-cell illustration have
  slightly different discretization results.

## New/updated files

- `numerics/value_of_information.py` — corrected, read-only VOI calculation
  against saved checkpointed known-S and conditional-red policies.
- `data/voi_grid.npz` — real known-time sample knots/costs.
- `data/voi_regret.npz` — broad and late conditional regret curve data.
- `diagnostics/session3_voi.json` — full numerical results and independent
  quadrature residuals.
- `diagnostics/session3_known_resolution.json` — cost convergence spots.
- `diagnostics/session3_known_curve.png`, `session3_regret_broad.png`,
  `session3_regret_late.png` — **diagnostics only**, not publication figures.
- `numerics/plot_session3_diagnostics.py` — regenerates those plots.
- `tests/test_session3_voi.py` — support/integration/known-reference regression.
- `results.json` — superseded VOI values removed from active current fields;
  Session 4/5/6 items remain clearly pending.

## Reproduction

```bash
# Populate any missing known-time reference knots (each point checkpointed):
python numerics/known_grid_checkpoint.py --start 22 --end 60 --step 1 \
    --starts rolling_0.65 rolling_0.95
python numerics/known_grid_checkpoint.py --times 26.5 27.5 28.5 29.5 \
    --starts rolling_0.65 rolling_0.95
# Challenge cached solutions using opposite-direction continuation:
python numerics/known_grid_checkpoint.py --start 35 --end 46 --reverse \
    --starts rolling_0.65

# NO expensive NLP in the following commands:
python numerics/value_of_information.py
python numerics/plot_session3_diagnostics.py
python tests/test_session3_voi.py
```

The per-start `max_cpu_time` in IPOPT is only a CPU quota; it is **not** a
reliable strict wall-clock watchdog for CasADi graph building or difficult
solve basins. Run sweeps in short separate subranges, using existing NPZ
checkpoints. Avoid all-in-one reoptimization over 39 times.

## Next sessions

Session 4: isolate acceleration *before the earliest possible green*, correct
its wording, assess realistic acceleration capping, and recompute any results
whose physical model changes. Session 5: publication figure redesign,
restructuring of numerical storyline, and possible fixed-mean width checks.
Session 6: article/boundary conditions/reproducibility audit and PDF rebuild.
