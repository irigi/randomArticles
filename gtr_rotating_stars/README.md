# Slowly Rotating Perfect-Fluid Spheres with an Adiabatic Gas Surface

Updated after the third review (revision 4). Start with the journal PDF below. The prior three-stage work is retained as historical calculation provenance.

## Article and supplement

- `output/pdf/final/rotating_star_journal.pdf`: 15-page RevTeX 4.2 article.
- `output/pdf/final/rotating_star_companion.pdf`: 14-page article using the requested IEEEtran/portable-fallback header. The fallback was used for this build.
- `output/pdf/final/rotating_star_supplement.pdf`: three-page benchmark and thermodynamic supplement.
- Matching `.tex` files are self-contained, apart from the adjacent `figures/` directory and `fontsetup.tex`.
- `article_body.tex` and `references.tex` are also supplied as editable components. `build_paper.py` assembles these components with `journal_preamble.tex` and `companion_preamble.tex` before compiling. Edit the components and run the build script to propagate changes to both layouts.
- `review3/CHANGES.md`: changes addressing the third review.
- `review3/AUTHOR_ACTIONS.md`: personal scientific review, disclosure details and public-deposit tasks.
- `review2/stage3/REVIEW_RESPONSE.md`: historical response to the second review.

The two article versions contain the same scientific text and figures. Equation numbering and layout differ. The author block is active. The requested companion layout is retained, with the font and language settings corrected.

## Build the PDFs

Install a normal TeX distribution with RevTeX 4.2, AMS packages, hyperref, PSNFSS `mathptmx`, English `babel`, URW Type 1 and RSFS fonts, graphicx, booktabs and the packages listed in each preamble. IEEEtran is optional for the companion; its portable fallback uses geometry and stfloats/placeins.

From the archive root:

```bash
python build_paper.py
```

This first assembles the shared article components and then runs three pdflatex passes for each document in its figure directory. No BibTeX run is required: references are embedded. The supplied PDFs are already built and visually checked. Different TeX distributions or the presence of IEEEtran may change page breaks.

## Reproduce calculations and figures

Use Python 3.12 and the pinned versions in `review2/stage2/requirements.txt`:

```bash
python -m pip install -r review2/stage2/requirements.txt
python review2/stage2/code/benchmark.py
python review2/stage2/code/thermodynamics.py
python review3/figures.py
```

These commands regenerate the new external benchmarks, thermodynamic data and seven updated figures. Rotation and radial-stability figures are retained from the earlier calculation. To regenerate those two as well, run `python regenerate_retained_figures.py`. Numerical scripts overwrite their generated data; make a copy before experimenting.

Root scripts implement the background, second-order rotation, fixed-EOS sequences, radial modes and reference-model comparisons. Their archived results are in `technical_review/` and `rewrite/`. The stage 1 audit, stage 2 report and stage 3 response are historical technical checkpoints. Their forward-looking to-do lists and old figure descriptions are superseded by review3/CHANGES.md.

All stellar calculations use G=c=1, with radius scaled to unity where stated. The benchmark polytrope has P=K epsilon^2. The equal-compactness comparison polytrope has P=K rho0^2 and epsilon=rho0+P. These are different EOS conventions. Constant-density benchmarks test I and Q only; their uncorrected finite-surface-density Love output is not used.

The published polytrope targets were extracted from the original vector curves of Yagi and Yunes, arXiv:1303.1528, Figs. 6 and 7 (right panels). The n=1 curve is Color11; legend text precedes its line sample. The archive retains three source vertices, coordinate calibration, source hashes and the published uniform-density table values. Full external papers and figures are not redistributed. To repeat extraction using full curves, obtain the authors' source from https://arxiv.org/src/1303.1528 and place IbarC-PRD.eps and QbarC-PRD.eps in review2/papers/yagi_source/. Without them the benchmark uses the retained vertices.

The supplied thermodynamic reference values are in `review2/attachments/thermodynamics/computed_values.json`. The main thermodynamic figure shows rest-mass density, enthalpy, internal energy and the effective polytropic coefficient. The supplement shows pressure per particle; it is not a core-temperature profile. An isentropic barotrope alone does not fix a unique finite-temperature EOS. The shallow outer Gamma minimum is real, followed by convergence to 5/3 from below.

## Validation and limitations

`QA_REPORT.json` records compilation and layout checks. `SHA256SUMS.txt` records all packaged file hashes. The supplement and `review2/stage2/data/` record external benchmark errors and graphical precision. The benchmark tests exterior quadrupole matching; empirical I-Love-Q fits remain comparisons, not correctness tests.

Admissibility boundaries and finite-parameter diagnostics remain numerical results, not interval-arithmetic certificates. The background is explicit in quadratures and frame dragging is algebraic. The second-order interior is represented by the regular Hartle system and a convergent integral series; a finite expression in standard special functions has not been obtained.

## Author actions before public release

Confirm author metadata and the AI-assistance disclosure, perform the personal scientific review needed to take responsibility for the claims, choose a target journal, and select code/data licensing terms. Consult review3/AUTHOR_ACTIONS.md, including model/version details and required funding/conflict statements. `zenodo_metadata_draft.json` and `CITATION.cff` are preparation files, not evidence of a deposit. No license grant, Zenodo deposit, DOI reservation, journal submission or arXiv posting has been made. After an actual deposit, add its DOI to the article and citation metadata. No repository DOI has been invented.

## Changes in this revision

Matese and Whitman (1980) is now cited and distinguished from the present dilute surface construction; the description of Whitman (1983) is grounded in its publisher abstract. Main displayed values are rounded without changing archived calculations. The PDFs use fully embedded fonts and English hyphenation. The disclosure explicitly covers substantive analytical, computational and writing assistance. No independent human audit or public DOI is claimed.
