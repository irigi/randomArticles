# Optimal Driving Toward a Red Traffic Light — Session 6 publication checkpoint

**Current published manuscript:** `article.tex` and compiled `article.pdf`. The four current vector figures are
`figures/fig1_known28_capped.pdf`, `fig2_information_capped.pdf`,
`fig3_anticipation_capped.pdf`, and `fig4_regret_capped.pdf`.
They are generated from **capped acceleration (4.136672736 m/s² net), 4 m/s² brake addition**, and the smooth `balance` resting convention.

**Current canonical numbers:** `results.json`. It has a single explicit publication model ID (`capped_b4_balance`).
The original uncapped and superseded calculations are retained in `data/`, diagnostics, and
`diagnostics/historical_results_through_session5b.json` for historical comparisons; they are not current-paper results. Legacy plot scripts and images reside only in `figures/legacy_uncapped/`.

## Paper build

```bash
pdflatex -interaction=nonstopmode -halt-on-error article.tex
pdflatex -interaction=nonstopmode -halt-on-error article.tex
```

## Figure reproduction from saved verified data

```bash
OPENBLAS_NUM_THREADS=1 python make_figures.py
```

## Verification (does not rerun all NLP optimizations)

```bash
OPENBLAS_NUM_THREADS=1 python tests/test_session3_voi.py
OPENBLAS_NUM_THREADS=1 python tests/test_session4_anticipation.py
OPENBLAS_NUM_THREADS=1 python tests/test_session5a_capped_voi.py
OPENBLAS_NUM_THREADS=1 python tests/test_session5b_publication_figures.py
OPENBLAS_NUM_THREADS=1 python tests/test_session6_article.py
```

For rebuilding computational results rather than only figures, see `SESSION1.md` through `SESSION5B.md`, especially the bounded independent commands in Session 5A. Older scripts under `numerics/` that summarize archived sessions should not be confused with generation of the current `results.json`.

## Semantics and caveats

- Brown curves = incomplete-information path **conditional on continued red**. Green and additional colors = post-green branches. Dashed charcoal = known-time reference when shown.
- `S=28` comparisons use the same physical realization; VOI integrates different realizations under each prior's own support.
- Anticipation counterfactual uses `u <= u_cruise(v)` **only for times strictly before 27 s**; the historical all-red restriction has a different numerical meaning.
- At rest, the solver balances modeled road load with tiny throttle and the fuel model incurs only idle fuel because `Gamma(0)=0`. This is a numerical convention, not literal advice to hold throttle at a stopped car.
- All reported optima are multistart **best-found feasible local candidates**; neither globality nor a general monotonic dependence on uncertainty width has been proved. The capped 0–100 calibration is not validation of a measured transient fuel map.

See `SESSION6.md` for paper revisions and verification results.
