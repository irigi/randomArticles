# Restricted driving policy: revised article

- `article.tex`: complete editable LaTeX source.
- `article.pdf`: compiled article.
- `figures/*.pdf`: vector figures used by the article.
- `figures/*.png`: convenient figure previews.
- `reproduce.py`: calculations and figure generation.
- `results.json`: numerical results from the supplied script.

Compile with `pdflatex article.tex` twice. The figures are already included.
To recompute, install NumPy, SciPy and Matplotlib, then run `python reproduce.py`.

The owner-measured fuel points are preserved exactly. Manufacturer specifications
refer to 2017 front-wheel-drive manual and DSG variants; the owner's exact
transmission and model year are unspecified. All remaining mechanical inputs and
the extension from cruise fuel use to transient fuel use are labeled assumptions.
The optimization is exact within the specified five-stage family, not claimed to
be the unrestricted global driving optimum. All monetary comparisons use the same
5 km route. The article contains source references and distinguishes homologation
fuel consumption from steady-speed measurements.
