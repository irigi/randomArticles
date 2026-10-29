# Slowly Rotating Perfect-Fluid Spheres with an Adiabatic Gas Surface

Revision 4. Start with `rotating_star_journal.pdf`.

## Layout

| Path | Contents |
|---|---|
| `rotating_star_journal.pdf` | 15-page RevTeX 4.2 article |
| `rotating_star_companion.pdf` | 14-page article in the IEEEtran/portable-fallback layout (fallback used for this build) |
| `rotating_star_supplement.pdf` | 3-page benchmark and thermodynamic supplement |
| `*.tex`, `figures/` | LaTeX sources and the figures they include |
| `build_paper.py` | Assembles the article sources and compiles all three PDFs |
| `code/` | Python scripts that reproduce every number and figure |
| `data/` | Results written by the scripts |
| `data/inputs/` | Fixed reference inputs that the scripts read but do not regenerate |
| `AUTHOR_ACTIONS.md` | Tasks for the author before submission |
| `CITATION.cff`, `zenodo_metadata_draft.json` | Draft release metadata (no deposit made yet) |

The journal and companion versions have the same scientific text and figures; equation numbering and layout differ.

## Build the PDFs

Install a TeX distribution with RevTeX 4.2, AMS packages, hyperref, PSNFSS `mathptmx`, English `babel`, URW Type 1 and RSFS fonts, graphicx, booktabs and the packages listed in each preamble. IEEEtran is optional for the companion; its portable fallback uses geometry and stfloats/placeins.

```bash
python build_paper.py
```

The shared text lives in `article_body.tex` and `references.tex`. The build script combines them with `journal_preamble.tex` and `companion_preamble.tex` into `rotating_star_journal.tex` and `rotating_star_companion.tex`, then runs three pdflatex passes for each document. Edit the components, not the assembled files. `--assemble-only` skips compilation. No BibTeX run is needed; references are embedded. Other TeX distributions, or having IEEEtran installed, may change page breaks.

## Reproduce calculations and figures

Use Python 3.12 with the pinned versions in `requirements.txt`:

```bash
python -m pip install -r requirements.txt
python code/run_all.py
```

This takes about a minute and overwrites `data/` and `figures/`; make a copy before experimenting. The scripts run in this order:

| Script | Produces |
|---|---|
| `technical_family.py` | Generating family and admissibility thresholds (`family_results.*`, `thresholds.json`) |
| `technical_stability.py` | Fixed-EOS sequences (`fixed_eos_*`) |
| `technical_radial.py` | Radial modes (`radial_*`) |
| `technical_verification.py` | Einstein-equation residuals |
| `rewrite_comparison.py` | Reference-model comparison |
| `benchmark.py` | External benchmarks (`benchmark_*`, `EPS_axis_calibration.json`, `published_vertices.json`) |
| `thermodynamics.py` | Thermodynamic profiles and checks, surface asymptotics |
| `figures.py` | Seven figures: background, family, thermodynamic profiles, universal relations, external benchmarks, gamma surface validation, pressure per particle |
| `figures_rotation_stability.py` | Rotation and stability figures, `surface_results.json` |

`check_star.py`, `drag_regular.py` and `drag_second.py` are shared modules for the background and second-order rotation. `finalize_audit.py` is a standalone check that prints an independent Volterra-series evaluation and matching checks. `data/rounded_values.tsv` maps computed values to the rounded values shown in the article.

## Conventions and data notes

All stellar calculations use G=c=1, with radius scaled to unity where stated. The benchmark polytrope has P=K epsilon^2. The equal-compactness comparison polytrope has P=K rho0^2 and epsilon=rho0+P. These are different EOS conventions. Constant-density benchmarks test I and Q only; their uncorrected finite-surface-density Love output is not used.

The published polytrope targets were extracted from the original vector curves of Yagi and Yunes, arXiv:1303.1528, Figs. 6 and 7 (right panels). The n=1 curve is Color11; legend text precedes its line sample. The project keeps three source vertices (`data/published_vertices.json`), the coordinate calibration and the published uniform-density table values. Full external papers and figures are not redistributed. To repeat the extraction from full curves, get the authors' source from https://arxiv.org/src/1303.1528 and place `IbarC-PRD.eps` and `QbarC-PRD.eps` in `data/inputs/yagi_source/`. Without them the benchmark uses the kept vertices.

The supplied thermodynamic reference values are in `data/inputs/thermodynamics_computed_values.json`. The main thermodynamic figure shows rest-mass density, enthalpy, internal energy and the effective polytropic coefficient. The supplement shows pressure per particle; it is not a core-temperature profile. An isentropic barotrope alone does not fix a unique finite-temperature EOS. The shallow outer Gamma minimum is real, followed by convergence to 5/3 from below.

## Validation and limitations

The supplement and the benchmark files in `data/` record external benchmark errors and graphical precision. The benchmark tests exterior quadrupole matching; empirical I-Love-Q fits remain comparisons, not correctness tests.

Admissibility boundaries and finite-parameter diagnostics are numerical results, not interval-arithmetic certificates. The background is explicit in quadratures and frame dragging is algebraic. The second-order interior is represented by the regular Hartle system and a convergent integral series; a finite expression in standard special functions has not been obtained.

## Before public release

Confirm author metadata and the AI-assistance disclosure, do the personal scientific review needed to take responsibility for the claims, choose a target journal, and select code/data licensing terms. See `AUTHOR_ACTIONS.md`. `zenodo_metadata_draft.json` and `CITATION.cff` are preparation files, not evidence of a deposit. No license grant, Zenodo deposit, DOI reservation, journal submission or arXiv posting has been made. After a deposit, add its DOI to the article and citation metadata.

## History

Earlier review rounds (audits, review responses, superseded scripts and diagnostics) were removed from the working tree to keep the project readable. They remain in git history before the commit that introduced this layout.
