# Checkpoint 4 — onset, clean counterfactual and mechanical realism

**Fourth of six checkpoints.** This is a **numerical / diagnostic** archive, not
an updated publication. The TeX body and current `article.pdf`, as well as the
four saved publication figures, have **not** been rebuilt; they still include
obsolete statements such as "full throttle starts at 26.625 s" and old VOI
values. A tiny figure *generator text* correction replaces "full throttle
begins" with "acceleration begins", but the existing PNGs were intentionally
not mixed with the new plots. Publication redesign and rebuild remain for
Sessions 5–6. Use this report/results.json for updated numerical values.

## Why two distinct acceleration constraints matter

The **legacy** counterfactual restricts the conditional-red path *throughout*
`[0,27.5]` to `u[k] <= u_cruise(v[k])`, preventing additional tractive force
above the cell-start steady-state requirement. It tests the economic benefit of
accelerating **anywhere on the still-red branch**. However, green first becomes
possible at 27 s, so it cannot by itself isolate *anticipatory acceleration
before any green is possible*.

The new independent `no_accel_before=27.0` restriction applies exactly the same
cruise-throttle ceiling only to cells fully before **27 s**, and leaves the
`[27,27.5]` portion free. `solve_uniform` requires that the cutoff be a node,
otherwise it throws `ValueError`; this avoids inadvertently restricting
post-cutoff cells. The constraint checks are repeated independently on the
saved solution, and `no_speed_increasing_throttle` cannot be combined with the
new option. Both constraints are *cellwise* direct-transcription approximations
and not global-optimality certificates.

## Updated expected costs, narrow `U(27,27.5)` prior (CZK)

| Vehicle model | Time value | Cells | Free acceleration | No increase until 27 s | No increase during all red | Specific pre-27 benefit | Entire-red benefit |
|---|---:|---:|---:|---:|---:|---:|---:|
| Original (uncapped, brake=7) | 300 CZK/h | 1100 | 3.613660859 | 3.613946101 | 3.614650714 | **0.000285242** | 0.000989855 |
| Original (uncapped, brake=7) | 900 CZK/h | 1100 | 9.398187546 | 9.402874757 | 9.408597243 | **0.004687211** | 0.010409697 |
| Capped (net accel 4.13667, brake=4) | 300 CZK/h | 1100 | 3.615397944 | 3.616782489 | — | **0.001384544** | — |
| Capped (net accel 4.13667, brake=4) | 900 CZK/h | 1100 | 9.400222554 | 9.410469040 | 9.419744607 | **0.010246486** | 0.019522053 |

Here "entire-red benefit" is the original, broader claim's numerical
counterfactual. It must **not** be presented as the value of pre-earliest-green
acceleration. The specific pre-27 benefit is the `before27 - unrestricted`
difference. The extra difference (`all_red - before27`) reflects permitting
acceleration during 27–27.5 s while the signal remains red, after green has
become possible.

At original mechanics and 900 CZK/h, the first cell with throttle above cruise
starts at **26.625 s** with `u=0.67627`, *not full throttle*. The next
cell at **26.650 s** has near-full throttle. At the capped mechanics and
1100 cells, near-full throttle begins around **26.425 s**. At 550 cells, the
throttle starts increasing at **26.40 s** and becomes near-full around
**26.45 s**. These small cell-level onset differences are discretization
artifacts: final prose should not imply sub-cell precision.

## Model realism and sensitivity separation

The original `A_POWER/(v+VREG)-r(v)` produces a net full-throttle acceleration
of **65.88228 m/s² at rest** (~6.72 g), with **6.75118 s** for 0–100 km/h.
This is an unusually aggressive low-speed model, not merely a somewhat optimistic
0–100 time. The capped model limits full-throttle net acceleration to
**4.136672736 m/s²** and gives **8.40000 s** for 0–100 km/h in this model.
It still represents a simplified acceleration envelope, **not** measured
engine transients or validated throttle-specific fuel consumption.

To disentangle cap from braking, the 550-cell **cap-only**, original `brake=7`
model (900 CZK/h) gives:

- unrestricted 9.400121150 CZK;
- no acceleration before 27 s 9.410359844 CZK;
- specific pre-27 benefit **0.010238694 CZK**.

The same comparison with cap **and** brake=4 gives unrestricted
9.400238901 and constrained 9.410476996, benefit **0.010238095 CZK** at 550
cells. Thus the cap, not the changed braking limit, dominates the altered
numerical advantage *within these tested solutions*.

The 550 → 1100 checks show the constrained 900-CZK/h cost moves by
`-0.0000262074` CZK for original mechanics and `-0.0000079564` CZK for capped;
unrestricted cap/900 shifts by about `-0.00001635` CZK. The main qualitative
pre-earliest-green result persists at both resolutions. As in Sessions 1–3,
IPOPT solutions are **local candidates**, not globally proved.

## Publication/model selection gate

**Recommendation:** the capped model deserves the *main* numerical treatment
on physical grounds, with the original as a historical comparison. However,
**do not silently switch** the main model yet. Session 3's known-time curve,
regret curve, and VOI apply to the **original uncapped model only**. Promoting
capped mechanics requires new capped broad/late unrestricted validation, a
new capped known-S grid and VOI calculation, and new figures/text. Given the
previous timeouts, this is better handled in **Session 5A** as a separately
resumable numerical checkpoint, then **5B** for graphics, and **6** for the
manuscript and PDF. If we instead keep the original model as main, make the
low-speed mechanical limitation far more prominent and report the capped
experiment as the robust/realistic comparison.

## Updated source, new data, and diagnostics

- `numerics/solve_policies.py` introduces `no_accel_before` with cutoff
  alignment validation and independent feasibility checks.
- `numerics/run_session4_case.py` solves **one** named variant + restriction +
  grid in each call, and writes `data/s4_*.npz` plus per-candidate JSON.
- `numerics/summarize_session4.py` rebuilds all values from saved NPZ and
  updates `results.json`/`diagnostics/session4_summary.json` without IPOPT.
- `numerics/plot_session4_diagnostics.py` generates the three diagnostic
  time-domain comparisons in `diagnostics/session4_*_comparison.png`.
- `tests/test_session4_anticipation.py` is a fast **no-NLP** validation of
  new constraints, cost ordering, cell onsets, calibration, cutoff guards and
  550/1100-cell consistency.
- `data/narrow900_noacc_N1100.npz` now uses the standardized `cost`/`lo`/`hi`
  schema, rather than the obsolete `J`/`low`/`high` schema. The old broad red
  constraint remains explicitly the **all-red** counterfactual, for backward
  reproduction; `make_figures.py` remains compatible with both cost keys.
- The PDF/manuscript main values, old mechanical-sensitivity policy arrays,
  known-S reference/VOI curves and official publication figures have not been
  changed. No stale value has been quietly reused for a different model.

## Short reproduction commands

```bash
# Independent bounded experiments; each overwrites only its own case files:
OPENBLAS_NUM_THREADS=1 python numerics/run_session4_case.py --case original900 --restriction before27 --cells 1100
OPENBLAS_NUM_THREADS=1 python numerics/run_session4_case.py --case capped900 --restriction before27 --cells 1100
# Fast, no IPOPT:
python numerics/summarize_session4.py
python tests/test_session4_anticipation.py
python numerics/plot_session4_diagnostics.py
```

For every publication comparison, verify `--case`, `--cells`, `--restriction`
and the recorded candidate diagnostics; do not compare numbers across
mechanical models as though they used the same objective/dynamics.
