# Checkpoint 2 — Uniform-prior / conditional-red solver robustness

This is the **second of six proposed sessions**. It builds on `SESSION1.md`.
The known-switch solver repairs from Session 1 remain intact and pass their
480-cell regressions after the changes in this checkpoint.

This checkpoint replaces the **unrestricted** saved policies for `U(0,60)`,
`U(26,30)`, and `U(27,27.5)` with results from a multistart uniform-prior NLP.
The manuscript, publication figures, VOI table/script output, mechanical
sensitivity conclusions, and the red-acceleration counterfactual are **not yet
revised** and must not be treated as re-certified by this archive.

## Solver repairs

`numerics/solve_policies.py` now gives `solve_uniform` the same safeguards as
the repaired deterministic solver:

- physically distinct starts: legacy, two rolling-delay starts, stop/wait,
  coast-then-stop, and brake/recovery;
- optional warm starts / horizon-rescaled continuation from other prior
  solutions;
- nonnegative velocity at the RK4 intermediate stages as well as at nodes;
- independent post-solve checks of dynamics defects and path bounds;
- per-start diagnostics and explicit selection of the lowest-cost feasible
  local candidate;
- a bounded `max_cpu_time` option and start-label subsets, because a poor
  high-dimensional start can otherwise consume an entire batch;
- `solve_uniform_family` for forward/reverse continuation studies over prior
  families.

The new `numerics/generate_session2_data.py` is deliberately resumable: each
paper-resolution policy can be regenerated and saved independently with
`--case broad`, `late`, `narrow300`, or `narrow900`.

## Publication-resolution unrestricted policies

Baseline mechanical model, smooth `rest_mode='balance'` convention:

| Prior / time value | Cells | Expected cost (CZK) | Minimum speed | Qualitative behavior |
|---|---:|---:|---:|---|
| `U(0,60)`, 300 CZK/h | 480 | 4.4078327398 | ~0 km/h | stops at line if red persists |
| `U(26,30)`, 300 CZK/h | 360 | 3.8430131661 | 28.346 km/h | rolling throughout red branch |
| `U(27,27.5)`, 300 CZK/h | 1100 | 3.6136608590 | 32.243 km/h | rolling; later recovery |
| `U(27,27.5)`, 900 CZK/h | 1100 | 9.3981875457 | 32.240 km/h | rolling; anticipatory recovery |

For the broad prior, five successful N=480 starts agree within roughly
`1.3e-7 CZK`; one deliberately slow rolling start can be difficult at this
resolution and is omitted from the bounded checkpoint generator after being
exercised at coarse resolution. For the late prior, all six N=360 starts
converge to the same solution to numerical precision. For each N=1100 narrow
case, the legacy and independent brake/recovery starts converge to the same
solution. Coarse six-start tests likewise find one basin for all three priors.
These are robustness checks, **not global-optimality certificates**.

Compared with the previous archive, the late and narrow trajectories are
unchanged to solver tolerance. The broad expected cost decreases by only
`2.3904e-6 CZK`; its maximum state changes are below 0.3 mm in distance and
0.0007 km/h in speed. Thus the uncertain-policy qualitative results survive
this repair even though the deterministic known-S results from Session 1 did
not.

## What the trajectories actually support

The robust policies confirm the reviewer's descriptive correction. At 10 s and
20 s, the broad-prior conditional-red trajectory is faster than the late-prior
one (about 44.09 vs 38.51 km/h and 38.57 vs 33.30 km/h). The distinction is
therefore **when delay is created**, not a demonstrated rule that "broader
uncertainty reduces the value of preserving speed." With `U(0,60)`, early green
is possible, so the policy initially preserves more speed and only commits to
full braking at about 24.625 s if red persists, reaching essentially zero speed
at about 26.0 s. With `U(26,30)`, substantial delay is guaranteed, so more of
that delay is created earlier while retaining a rolling approach through the
late window.

The two priors change mean, width, earliest green, and latest green at once.
They therefore do **not** isolate an uncertainty-width effect. Any such claim
must wait for a fixed-mean width sweep (scheduled for the later numerical/
figure session if still useful).

## Rest convention

The publication solver keeps the smooth legacy force-balance convention for
this checkpoint. Literal `u=b=0` at `v=0` would give negative acceleration in
that ODE. The optional unilateral `rest_mode='contact'` introduced in Session 1
is physically clearer but nonsmooth at the boundary and was more
initialization-sensitive in the broad-prior NLP.

Importantly, under this paper's calibrated fuel law `gamma(0)=0`. The small
balancing throttle used by the smooth ODE at exactly zero speed therefore adds
no traction-fuel term above the idle fuel term. It is a numerical convention
for holding the state, not a claim about literal pedal position, and should be
said explicitly in the final manuscript.

## Commands

```bash
# Fast/coarse multistart regression
python numerics/check_uniform_multistart.py

# Paper-resolution unrestricted policies, independently resumable
python numerics/generate_session2_data.py --case broad
python numerics/generate_session2_data.py --case late
python numerics/generate_session2_data.py --case narrow300
python numerics/generate_session2_data.py --case narrow900

# Recheck Session 1 after the solver edits
python numerics/check_known_multistart.py
```

Detailed candidate outcomes are in
`diagnostics/session2_uncertain_solver.json`.

## Still pending

1. Recompute the dense known-S curve, correct the VOI integration script, and
   build the regret curve (Session 3).
2. Replace the current anticipatory-acceleration counterfactual with one that
   isolates acceleration before the earliest possible green; revisit the
   mechanical model/sensitivity and exact acceleration-onset wording (Session 4).
3. Rebuild publication figures and restructure the numerical storyline only
   after those numbers stabilize (Session 5).
4. Rewrite the article, boundary-condition discussion, and final conclusions
   from the verified results (Session 6).
