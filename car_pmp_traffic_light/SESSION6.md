# Session 6 — Final article revision and validation

The **capped** numerical model (`capped_b4_balance`, net full-throttle acceleration cap 4.136672736 m/s², brake addition 4 m/s², balance-at-rest convention) is now the **principal paper model**. Legacy original/uncapped comparisons are explicitly historical and kept separately. The seven-page `article.pdf` was built twice from `article.tex`, and four figures are included as **vector PDFs** generated in Session 5B.

## Scientific corrections

- The paper's central question is how to allocate sufficient delay while preserving speed and position at green. It begins with the cruising and free-coasting benchmarks (21.6 and 25.30 s).
- Mechanical parameters and fuel-law re-normalization are specified for the capped case, with the uncapped 65.88 m/s² low-speed artifact candidly noted.
- Green terminal value, crossing vs arrival, finite-fine crossing boundary `K+V_G`, prohibited-crossing viability, waiting-at-rest balance convention, time-support terminal link, red-history contingent policy, HJB control-switching and survival objective are explained.
- Known-S brake–coast–recover example leads the figures. Broad and late priors are interpreted in terms of **when speed is sacrificed**, not as a controlled uncertainty-width experiment.
- The narrow-prior counterfactual forbids speed-increasing throttle **only before 27 s**, rather than on the whole red branch; its cost benefit is 0.010246 CZK at 900 CZK/h in the capped model.
- Information values are computed consistently under the capped model and plotted as prior-support-specific regret; old uncapped values are mechanical sensitivity results only.
- The deterministic solver's demonstrated local-minimum failure and multistart repair are stated, with explicit numerical limitations rather than claims of global optimality.

## Principal values (CZK)

| Metric | Broad U(0,60) | Late U(26,30) |
|---|---:|---:|
| Expected capped feedback cost | 4.463750626 | 3.850498274 |
| Expected capped perfect-info reference | 4.278767726 | 3.673527570 |
| Estimated capped VOI | 0.184982900 | 0.176970704 |
| S=28 realization feedback cost | 4.133089972 | 3.849904126 |

The capped known-S=28 reference cost is 3.675504388 CZK. Capped narrow-prior 900 CZK/h unrestricted / no-pre27 costs are 9.400222554 / 9.410469040 CZK, respectively. These are best-found feasible local numerical candidates.

## Verification

- Article compiled with two `pdflatex` passes; **7 A4 IEEE double-column pages**.
- New figures inspected on pages 4–7 after rendering PDF to images; caption numbers and references verified.
- All numerical-regression checks run: `test_session3_voi.py`, `test_session4_anticipation.py`, `test_session5a_capped_voi.py`, `test_session5b_publication_figures.py`, `test_session6_article.py`.
- `results.json` now clearly identifies the capped publication model, while superseded original/checkpoint contents are archived under `diagnostics/historical_results_through_session5b.json` and the prior `SESSION*.md` files.
- Old unnumbered/uncapped figure compatibility copies removed from `figures/` top level; originals preserved under `figures/legacy_uncapped/`.

Remaining limitations: no global optimality proof; no systematic fixed-mean prior-width sweep, no transient fuel-law calibration, finite-fine numerical exploration or comfort modeling. To regenerate solver outputs, consult bounded checkpoint scripts (not all commands are fast). The publication values refer to the supplied saved solutions.
