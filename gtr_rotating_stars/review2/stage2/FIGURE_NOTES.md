# Stage 3 figure integration notes

All seven figures are supplied as vector PDF and 180-dpi PNG. Main panels are
designed at about 181 mm width. Final numbering belongs to Stage 3.

| File | Placement | Required caption information |
|---|---|---|
| `background.pdf` | Replace Figure 1 | Density and pressure, reconstructed EOS, sound speed, and adiabatic index for alpha=0.1 and 0.3. The inset resolves the outer local minimum and approach to Gamma=5/3 from below. The endpoint is an analytic limit; interior points use the hydrostatic pressure integral near the surface. The reference powers are normalized comparisons, not global EOS fits. |
| `thermodynamic_profiles.pdf` | New main figure after EOS reconstruction | Four panels: ideal-gas-equivalent temperature, enthalpy per particle, internal energy above rest mass, and effective polytropic coefficient. Define rho0=m_b n, H_b=(epsilon+P)/n, and K_eff=P/rho0^(5/3). T_ideal=P/(n k_B) is a conditional interpretation. The barotrope alone does not determine physical temperature or heat capacity. |
| `universal.pdf` | Replace Figure 4 | I-Q relation and fractional residuals in I relative to the Yagi-Yunes I-Q and I-Love fits. Add Tolman VII, rest-mass n=1 polytrope, and Buchdahl points at the two common compactnesses. Family points cover 0.1<=alpha<=8; connecting segments guide the eye, not continuous fitted results. Hollow blue family circles denote the sampled radially unstable models (alpha=1,3,8). Reference-model symbols are hollow solely for legibility. These fits are comparisons, not correctness tests. |
| `family.pdf` | Replace existing family plot | Vertical lines mark monotone-stiffening, radial-mode and causality thresholds. Thresholds come from root-finding and radial maximization, not interpolation of plotted markers. The largest-alpha points are neither all radially stable nor all monotonically stiffening. Numerical boundaries are not interval certificates. |
| `external_benchmarks.pdf` | Supplement or compact validation appendix | Left: differences from Chandrasekhar-Miller (1974), Table I, for uniform-density stars. Right: differences from digitized n=1 curves in Yagi-Yunes (2013), Figs. 6 and 7. Error bars show a conservative graphical-resolution allowance including both coordinate directions, not statistical errors or published solver uncertainties. The n=1 benchmark uses P=K epsilon^2. |
| `gamma_surface_validation.pdf` | Supplement | The shallow outer dip is real. The linear asymptote is Gamma=5/3-5z/9+O(z^2), z=sqrt(1-r^2/R^2). For alpha=0.3 this outer minimum is local; the centre has a lower Gamma. |
| `temperature_profile.pdf` | Supplement; redundant in main text | Normalized and dimensionless ideal-gas-equivalent temperature. Do not infer a Kelvin scale without specifying particle mass and the extra thermal interpretation. Mechanical equilibrium does not impose Tolman thermal equilibrium here. |

The supplied Figure 1 and thermodynamics archives have been incorporated by
recomputing the profiles with consistent central values, more stable surface
derivatives, and shared styling. No supplied curve was smoothed or forced flat.
The main article and its figure numbering have not yet been changed.
