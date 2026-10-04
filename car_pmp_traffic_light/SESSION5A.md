# Session 5A — capped mechanics, known-time reference and value of information

This is a **numerical checkpoint only**. The publication `article.tex`, `article.pdf`
and four figures remain from the previous, **uncapped** main example. Do not cite their
VOI table or captions as capped-model results. Session 5B will regenerate plots;
Session 6 will revise the manuscript and rebuild the PDF.

## Selected model

The same calibrated **net full-throttle acceleration cap** of 4.136672736 m/s²
from Session 4, 4 m/s² maximum brake parameter, 300 CZK/h time value,
and smooth `rest_mode=balance` convention (not literal pedal behavior).
`Model.H(v)` is recomputed for this dynamics model. No uncapped green cost or
known-time references are reused.

## Capped model (CZK, best feasible *local* candidates)

| Prior | Feedback expected | Perfect info expected | VOI | Independently integrated regret |
|---|---:|---:|---:|---:|
| Broad U(0,60) | 4.463750626 | 4.278767726 | **0.184982900** | 0.184980174 |
| Late U(26,30) | 3.850498274 | 3.673527570 | **0.176970704** | 0.176970678 |

Both priors' feedback objective is independently integrated by cellwise Gauss
quadrature; differences from NLP are **-2.73e-06** CZK
and **-2.58e-08** CZK, respectively.

## Original model vs capped mechanics

| Prior | Original VOI | Capped VOI |
|---|---:|---:|
| U(0,60) | 0.132081 | 0.184983 |
| U(26,30) | 0.171021 | 0.176971 |

The ordering changes: in the capped model, **broad > late**, while with the
original uncapped mechanics, **late > broad**. These two priors also differ
in mean and support, so this is **not** evidence for a universal trend with
uncertainty width. The broad-policy cost is more sensitive to recovery from
stopping when low-speed acceleration is capped.

## Reference, robustness and provenance

- 43 independently checkpointed capped known-`S` candidates, integers 22..60
  and half-second values 26.5,27.5,28.5,29.5 at **N=480**; plus the constant
  free-flow optimum on `[0,21.6]` (66 interpolation knots in total).
- Forward continuation plus an independently challenged **18-point** reverse
  subset, using physically interpretable rolling starts and cached solutions.
- Capped known costs at S=28,40,60 s: **3.675504388**,
  **4.960873772**, **6.835906248** CZK.
- At S=60 s, N=480 yields 6.835906248 CZK and N=960 yields
  6.835870706 CZK, a difference of
  **-0.0000355 CZK**.
  This is a spot check, not a dense convergence or globality certificate.
- `data/u_0_60_capped.npz`, `data/u_26_30_capped.npz` and their solver diagnostics
  are explicitly model tagged. Existing uncapped files, original-session VOI
  plots and tables are preserved unchanged for historical reproduction.
- `data/voi_capped_grid.npz`, `data/voi_capped_regret.npz` and
  `diagnostics/session5a_voi.json` are model-specific. No stale references used.

### Reproduction

```bash
OPENBLAS_NUM_THREADS=1 python numerics/session5a_uniform.py --case broad
OPENBLAS_NUM_THREADS=1 python numerics/session5a_uniform.py --case late
OPENBLAS_NUM_THREADS=1 python numerics/session5a_known.py --direction forward
OPENBLAS_NUM_THREADS=1 python numerics/session5a_known.py --direction reverse --times 22 25 26 26.5 27 27.5 28 28.5 29 29.5 30 33 35 40 45 50 55 60
OPENBLAS_NUM_THREADS=1 python numerics/session5a_known.py --times 60 --cells 960
python numerics/session5a_voi.py
python numerics/plot_session5a.py
python tests/test_session5a_capped_voi.py
```

**Limitations:** these costs are best found IPOPT *local* candidates, not
certified optima. The 480-cell reference, interpolation between switch times,
and stage dynamics introduce numerical error, though independent quadrature,
multistart, reverse challenges and the 960-cell endpoint test bound selected
sources. No fixed-mean uncertainty-width sweep has been performed.
