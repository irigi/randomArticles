# Article revision - 4 October 2026

This revision changes the manuscript and its compiled PDF only. All original numerical code, tests, data, figure files, historical diagnostics, and other source-package files are retained byte-for-byte. This note is the only added file in the source archive. Captions and figure placement are part of the manuscript and have been edited; the image and vector-figure files themselves have not been changed.

## What changed in the article

- Rebuilt the explanation around creating delay while preserving position and speed for the remaining journey.
- Put route constraints, vehicle dynamics, fuel use, and the objective in one model section.
- Completed the green-continuation verification and derived an equivalent objective that expresses the additional cost of delay, with survival weighting under uncertainty.
- Explained the shared red path before introducing the HJB and control-switching conditions.
- Reorganized the numerical narrative around four questions: known-time delay allocation, uncertain-time slowing, acceleration before the earliest green, and realized information regret.
- Added a small table of state at green, destination arrival time, and fuel use for the common realization S = 28 s. These values were derived from the existing saved trajectories during the review; no optimizations were rerun.
- Explained the peak and subsequent decline in broad-prior regret.
- Corrected the description of stationary control degeneracy, the scope of numerical checks, and the confounding of acceleration and braking in the alternative-model comparison.
- Moved finite red-crossing penalties to a short appendix and consolidated numerical qualifications.

## Follow-up changes outside the manuscript

### Needed if the model is to represent one specific Octavia variant

The existing results combine a drag coefficient of 0.302, listed for the Combi, with an 8.4 s 0-100 km/h benchmark for the liftback. The cited liftback has drag coefficient 0.294; the Combi acceleration time is 8.5 s. The revised text openly treats the existing parameters as an illustrative set. To make a specification-consistent model, choose one variant, update the corresponding mechanics and acceleration calibration, and regenerate the policies, numerical summaries, and figures. None of those changes has been made here.

### Needed for stronger causal or numerical claims

- To isolate the effect of the acceleration cap on VOI, hold brake strength fixed. The current main and uncapped comparison models use brake additions of 4 and 7 m/s^2, respectively. The narrow-prior archive contains a limited cap-only comparison, but the VOI comparison is confounded.
- To claim an effect of uncertainty width, add a controlled prior family, such as a fixed-mean sweep. The two existing priors illustrate different timing opportunities.
- To prescribe zero brake while waiting, impose a canonical stationary control convention in the solver. At present it permits simultaneous throttle and brake satisfying force balance. The rewritten text describes this accurately; the idle-fuel cost at exact rest is unaffected.
- Global optimality and quantitative error bounds would require additional analysis. The current checks support the reported candidates and their evaluation, not certified optima.

### Optional figure improvements

- Figure 2: a speed-versus-time panel would directly support the discussion of when speed is sacrificed. The current image remains usable with the revised caption and added table.
- Figure 4: annotate the feedback braking/stopping transition and the regret peak near 26.5 s. Its existing curves already support the revised interpretation.
- Round monetary annotations in a future figure refresh to match the precision used in the prose. The detailed values in the unchanged images remain consistent with the saved calculations.

### Documentation and tests

The existing README, session notes, and tests are preserved. Their historical claims, section numbering, exact manuscript strings, and layout assumptions may no longer describe the rewritten article. In particular, manuscript-specific assertions in `tests/test_session6_article.py` should be reviewed before treating that test as a gate for this revision. No code or tests were edited or rerun for the rewrite. Build the revised `article.tex` with the existing class and figure files using two pdflatex passes.

## Verification used in the article

The preceding review independently reintegrated the saved piecewise-constant controls using SciPy DOP853, respecting control-cell boundaries, with relative tolerance 2e-11 and absolute tolerance 1e-12. It used the existing smooth force-balance dynamics without clipping state trajectories. The cost differences (reintegrated minus saved, CZK) were:

| Saved policy | Difference |
| --- | ---: |
| Known 28 s, 480 cells | -0.000003207 |
| Broad U(0,60), 480 cells | -0.000004913 |
| Late U(26,30), 360 cells | -0.000000593 |
| Narrow, 900 CZK/h, unrestricted, 1,100 cells | -0.000000850 |
| Narrow, 900 CZK/h, no pre-27 acceleration, 1,100 cells | -0.000001532 |
| Narrow, 300 CZK/h, unrestricted, 1,100 cells | -0.000002647 |
| Narrow, 300 CZK/h, no pre-27 acceleration, 1,100 cells | -0.000004589 |

These are checks of fixed saved policies, not new optimization results. The model's continuation function was also compared with direct numerical quadrature during the review. The manuscript rewrite was compiled, checked for unresolved references and overfull boxes, and visually inspected. Package verification compares every original archive member: only `article.tex` and `article.pdf` are replaced.
