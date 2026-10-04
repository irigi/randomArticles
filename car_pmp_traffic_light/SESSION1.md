# Checkpoint 1 — Known-switch deterministic solver and regression checks

This is the **first of six proposed sessions**, not a revised publication. The
article PDF, `article.tex`, figures, prior data, `results.json`, and VOI script
are carried forward unchanged and are **not certified** by this checkpoint.
Do not cite the old VOI table as recalculated. Uncertain-policy solver, VOI,
main figures and narrative are scheduled for later sessions.

## Repairs

- `numerics/solve_policies.py` now runs physical initializations (legacy,
  rolling at two speed targets, stop/wait, and coast-then-stop), records
  individual successes/failures and their costs, and selects the cheapest
  independently feasible local candidate. Warm starts from adjacent switch
  times permit **both forward and reverse continuation**. Start subset and
  per-start IPOPT CPU budgets can be specified to keep experiments bounded.
- Nonnegative velocity is enforced both at transcription nodes and at three
  intermediate RK4 velocity stages. Since `d` decreases with nonnegative
  velocity, `d>=0` at nodes prevents the numeric red branch from crossing the
  line before the switch. Path residuals and speed bounds are checked again
  numerically after optimization.
- `numerics/model.py` now explicitly labels its legacy **force-balance at
  rest** (`rest_mode='balance'`, default, for reproducibility). With `v=0`,
  this model's raw moving ODE needs balancing traction and/or braking to hold
  zero speed. An optional `rest_mode='contact'` gives a unilateral reaction
  at rest, so `u=b=0` holds position without fuel-consuming traction. The
  rolling S=40 and S=50 optimum costs remain the same in that variant, as
  does the S=60 stop candidate. This is not yet a complete uncertainty-model
  revalidation; Session 2 should select the common model for all branches.
- `save_npz` retains only numeric policy keys and simple text metadata,
  preventing `candidates` diagnostic dictionaries from leaking into the
  archival `.npz` data as pickle-dependent object arrays.

## Direct 480-cell regression points, baseline parameters (CZK)

| Known green S | Former single start | Multistart | Minimum speed |
|---:|---:|---:|---:|
| 40 s | 5.06083288 (stops) | 4.95568147 | 17.706 km/h |
| 50 s | 5.94712354 (stops) | 5.91244709 | 10.347 km/h |

The S=40 cost agrees to the shown decimals with the independent review.
The S=50 cost differs by approximately 0.00006 CZK from the reviewer's
5.91239 CZK value; this difference requires closer study of formulations,
quadrature, or local trajectories before claiming exact reproduction.

Coarse checkpointed **forward/reverse** sweep at N=480 (see
`diagnostics/known_sweep`, not a full 1-second VOI grid):

| S (s) | Lowest candidate cost (CZK) | Minimum speed (km/h) |
|---:|---:|---:|
|25|3.311423034|38.580|
|28|3.673977424|31.410|
|35|4.448381722|22.503|
|40|4.955681470|17.706|
|45|5.441164337|13.717|
|50|5.912447093|10.347|
|55|6.373883656|7.387|
|60|6.823812573|0|

The 60-second stop may be optimal or a local candidate; **this sweep does not
certify global optimality**. Nor do these findings prove all later switches
are rolling. What has been established is that the former blanket claim
that known late switching *forces* a full stop is false at S=40/50.

## Commands

```bash
python numerics/check_known_multistart.py
python numerics/known_grid_checkpoint.py --times 25 28 35 40 45 50 55 60
python numerics/known_grid_checkpoint.py --times 25 28 35 40 45 50 55 60 --reverse
# For subsequent incremental batches of the full 1-second known-switch grid:
python numerics/known_grid_checkpoint.py --start 22 --end 30 --step 1
```

The checkpoint sweep is resumable. It saves results as separate files after
**each** switch time instead of putting the entire campaign at risk of a
session timeout. Reverse mode challenges each solution again with alternate
warm starts and keeps the lower cost. For scientific comparison, keep model,
cell count, cost quadrature and prior fixed across the sweep.

## To be addressed later

1. Finish dense, N-converged known-S grid and corrected perfect-information
   integrals (Session 3). The present VOI script is known to contain a mask
   error and its sample-interval description is inconsistent.
2. Cross-check the rest/contact convention with the uncertain-policy NLP
   and green continuation consistently (Session 2).
3. Verify anticipatory counterfactual, mechanical model and main figures
   before rewriting the paper (Sessions 4–6).
