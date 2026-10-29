# Changes after the third review

## Article and interpretation

- Restored Jan Olšina's author name, affiliation, email and ORCID in the main article; restored the author and affiliation in the supplement.
- Added a data-and-code availability statement referring to the supplied reproducibility archive. No public repository or DOI is claimed.
- Added an extensive AI-assistance disclosure, a computational-provenance subsection, and attribution of AI-assisted plotting code in figure captions. The disclosure includes analytical development, derivations, code, literature work, interpretation and drafting. It does not claim that the author has already independently rederived every result.
- Replaced the main thermodynamic figure's temperature panel with normalized rest-mass density. The supplement retains the pressure-per-particle curves, explicitly labelled as such, without assigning a core temperature.
- Explained that the leading surface relations are compatible with a thermal isentropic monatomic ideal gas or a cold noninteracting degenerate Fermi gas. The stiff interior is outside those free-gas descriptions; interactions are a possible explanation, not an inferred microscopic EOS.
- Related the comparatively low radial-stability threshold to the decline in central adiabatic index across the generating family. The text distinguishes this interpretation from the eigenvalue calculation and from a fixed-EOS sequence.
- Rounded main-text and main-table computed values to at most six significant figures. Full numerical data are unchanged. The first-law residual is explicitly evaluated from unrounded data.

## Literature check

Matese and Whitman, *Physical Review D* **22**, 1270–1275 (1980), DOI [10.1103/PhysRevD.22.1270](https://doi.org/10.1103/PhysRevD.22.1270), was read from the attached six-page scan. Its actual title is *New method for extracting static equilibrium configurations in general relativity*. The gas-motivated example in Sec. IV is an approximation to the degenerate Fermi-gas EOS at high density, not a regular dilute-gas surface model. The paper identifies the example as singular, gives a finite surface density in Eq. (4.13), and states a nonzero surface sound speed squared of 2/9. The introduction now cites both its generating method and this physical distinction. The attached paper is not redistributed in the code archive.

Whitman (1983), DOI [10.1103/PhysRevD.27.1722](https://doi.org/10.1103/PhysRevD.27.1722), was checked against the publisher's abstract. The introduction now describes the regular-interior/generalization result at that documented level, separately from the first-order rotating solutions in the previously inspected 1985 paper. This revision does not claim that the full 1983 text was inspected.

## Typesetting

- Replaced `times` with `mathptmx`; added English `babel`.
- Added explicit pdfTeX font-map overrides for embedded URW Times-compatible outlines, together with required mathematical fonts. `pdffonts` confirms embedding for every font in all three output PDFs.
- Set English left/right hyphenation minima to 2/3. A TeX probe gives `ex-pres-sions`, `ex-act`, `bench-marks`, `in-de-pen-dent`, rather than the previously reported one-letter tails.
- Preserved the requested IEEEtran/fallback companion layout alongside the RevTeX journal version.
- Visually inspected the rebuilt article and supplement. Compilation has no unresolved references/citations or overfull boxes.

The changes clarify presentation and interpretation; they do not alter the EOS, metric construction, numerical datasets or external benchmark values.
