# Stage 2: validation and figure checkpoint

17 October 2026. Stage 2 of 3 is complete. This package contains the new calculations and figures; the article rewrite is reserved for Stage 3.

## 1. External quadrupole validation

The comparison solver now reproduces two independent published calculations. These tests use the matched exterior quadrupole, rather than an empirical universal-relation fit.

### Uniform-density table

Source: S. Chandrasekhar and J. C. Miller, *On slowly rotating homogeneous masses in general relativity*, MNRAS **167**, 63–80 (1974), [doi:10.1093/mnras/167.1.63](https://doi.org/10.1093/mnras/167.1.63), Table I, p. 73.

The exact interior Schwarzschild background is used at R=1. If C=M/R, then E=8 pi epsilon=6C, f=1-2Cr^2, b=sqrt(1-2C), a=sqrt(f), exp(Phi)=(3b-a)/2 and Pi=8 pi P=6C(a-b)/(3b-a). Only I and Q are retained; the comparison routine’s uncorrected Love output is discarded because this model has a surface density jump. No mass-correction benchmark is claimed.

For R_S=2M, the published inertia column is I/R_S^3 and its positive quadrupole column is Q_positive/(J^2/R_S). Hence Ibar=8 I_table and Qbar=Q_table/2. Our signed oblate quadrupole is negative, Qbar=-QM/J^2.

| C | Published Ibar | Calculated Ibar | Difference (%) | Published Qbar | Calculated Qbar | Difference (%) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.25000000 | 8.392000 | 8.39506760 | +0.03655 | 3.262500 | 3.26301980 | +0.01593 |
| 0.16666667 | 16.984000 | 16.98235350 | -0.00969 | 7.835000 | 7.83589182 | +0.01138 |
| 0.10000000 | 43.896000 | 43.89661657 | +0.00140 | 18.865000 | 18.86252015 | -0.01315 |

Maximum absolute differences are 0.0366% in Ibar and 0.01594% in Qbar. These are historical numerical targets: agreement is not uniformly within half the last printed digit, so the paper should quote the measured percentage agreement rather than claim exact reproduction of every rounded entry. Reducing the central cutoff s0 from 1e-7 to 1e-8 changes these results by less than 2e-13 relatively.

### Published index-one polytropic curves

Source: K. Yagi and N. Yunes, *I-Love-Q relations in neutron stars and their applications to astrophysics, gravitational waves, and fundamental physics*, PRD **88**, 023009 (2013), [doi:10.1103/PhysRevD.88.023009](https://doi.org/10.1103/PhysRevD.88.023009), [arXiv:1303.1528](https://arxiv.org/abs/1303.1528). The right-hand panels of Figs. 6 and 7 supply Ibar(C) and Qbar(C).

The benchmark uses their EOS P=K epsilon^2, not our table’s different relativistic rest-mass polytrope P=K rho0^2, epsilon=rho0+P. With scaled Pi=E^2, the logarithmic enthalpy is h=2 log(1+E), allowing a regular TOV integration. All three benchmark configurations have causal central sound speeds.

**Extraction correction to Stage 1:** the n=1 curve is Color11 (orange solid), not Color10. In these EPS files, the legend text precedes its corresponding line sample. Color10 is Shen. The initial mismatch was resolved by inspecting the rendered source figure; the wrong curve is not used in any final result.

The saved source vertices avoid interpolation of the published curve. Their EPS coordinates have four decimal places. The logarithmic world windows are C=0.004…0.6, Ibar=3…1000 and Qbar=0.7…300, with viewport x=0.15…1.15 and y=0.15…0.85. Ten checked tick coordinates agree within 4.80e-5. We allow +/-1e-4 per coordinate for rounding plus calibration and propagate both C and vertical-coordinate effects. This conservative graphical allowance is about 0.16–0.20%, not a rigorous interval or an uncertainty in the original numerical integration.

| C | Digitized Ibar | Calculated Ibar | Difference (%) | Digitized Qbar | Calculated Qbar | Difference (%) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.07756275 | 47.882190 | 47.88931381 | +0.01488 | 20.559615 | 20.55427473 | -0.02598 |
| 0.15409115 | 13.562786 | 13.56418359 | +0.01031 | 6.289881 | 6.28802543 | -0.02951 |
| 0.20005621 | 8.692976 | 8.68916859 | -0.04380 | 3.576748 | 3.57551423 | -0.03448 |

All six values agree within the graphical allowance. Maximum differences are 0.04380% in Ibar and 0.03449% in Qbar. Central-cutoff refinement changes the selected results by less than 2e-11 relatively. This is a comparison against published curves, not against the Yagi-Yunes fit.

Finally, applying the externally tested comparison routine to our own backgrounds reproduces the stored main-solver I coefficients within 5.1e-11 relatively and Q coefficients within 4.5e-10 at alpha=0.1 and 0.3. This last check is internal, but establishes the connection between the benchmarked routine and the manuscript results. None of these tests constitutes a global proof of the Einstein equations; the earlier local residual and convergence checks remain complementary.

## 2. Reference-model residuals

The revised universal-relation figure now includes all comparison models. The residual convention is Ibar/Ibar_fit-1 for both fits; Lambda=2 k2/(3 C^5). The rest-mass polytrope remains explicitly distinguished from the energy-density polytrope used above.

| Reference model | C | I-Q residual (%) | I-Love residual (%) |
|---|---:|---:|---:|
| Tolman VII | 0.08207711 | +1.5892 | +2.0784 |
| Rest-mass n=1 polytrope | 0.08207711 | -0.1410 | +1.3271 |
| Buchdahl gas | 0.08207711 | +1.2130 | +1.8901 |
| Tolman VII | 0.15884974 | -0.5387 | +0.5374 |
| Rest-mass n=1 polytrope | 0.15884974 | -2.3089 | -0.1770 |
| Buchdahl gas | 0.15884974 | -0.1033 | +0.6903 |

At C=0.15885 the family’s I-Q residual is -2.8258%, versus -2.3089% for the rest-mass n=1 reference. The near-zero residual at alpha=0.1 is not presented as validation.

## 3. Thermodynamics and the surface

The supplied supplement is reproduced, with a regular exact central limit and cancellation-resistant surface evaluation. Let H=E+Pi, z=sqrt(1-r^2/R^2), a_j=d(log j)/ds and g=Phi_r/r in R=1 units. The implementation uses

    H=-4 f a_j,  dPi/dz=z H g,
    h(z)=integral_0^z t g(t) dt,
    rho0_scaled=H exp(-h),  H_b/m_b=exp(h),
    theta=Pi/rho0_scaled,  u_int=expm1(h)-theta.

Differentiating H directly gives dE/dz=dH/dz-dPi/dz, avoiding the subtraction of tiny surface-density derivatives in the original background output. Pressure is obtained from the hydrostatic integral when z<0.2. expm1(h) avoids cancellation in the internal energy.

The independently checked surface expansions, with g_b=C/(1-2C), are

    Gamma=5/3-(5/9)z+O(z^2),
    theta=(g_b/5)z^2[1-z/3+O(z^2)],
    u_int/theta -> 3/2.

Thus the curves really dip below 5/3 before approaching it. At z=1e-4 the computed Gamma values are 1.666611149 and 1.666611151, while u_int/theta is 1.50008330 for both examples. The plotted endpoint is the analytic limit; the approach is resolved with computed interior points.

| alpha | theta_c | H_b,c/m_b | u_int,c | Outer local Gamma minimum | r/R at that minimum |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 0.0556364488 | 1.1081168593 | 0.0524804105 | 1.6411901770 | 0.9947475660 |
| 0.3 | 0.1514745706 | 1.3337462014 | 0.1822716308 | 1.6431289177 | 0.9956410412 |

For alpha=0.3 the outer minimum is not the global minimum: Gamma_c=1.53015. The finite-difference check of d(log rho0)/dr=(dE/dr)/(E+Pi) has scaled residual below 1.7e-10 at six interior radii. Algebraic pressure and an independently integrated hydrostatic pressure agree within 6.4e-10 relatively at the three checked inner z values. Metric and hydrostatic logarithmic enthalpies agree within 3e-13 absolutely at the checked points. These are numerical checks, not interval-certified bounds.

**Interpretation:** theta=P/(m_b n) equals k_B T_ideal/m_b only if the ideal-gas relation P=n k_B T is additionally adopted. The barotrope does not specify a unique physical temperature or heat capacity. H_b=(epsilon+P)/n is enthalpy per particle; it tends to m_b at the surface. Its rest-mass-subtracted value and logarithmic enthalpy vanish. At finite entropy it is not generally the Gibbs chemical potential. Mechanical equilibrium does not assert full Tolman thermal equilibrium.

## 4. Figure and manuscript handoff

Seven vector figures and PNG previews are ready. `FIGURE_NOTES.md` specifies placements and caption requirements. Figure 1’s label is visible and its surface inset resolves the limiting behavior. The thermodynamic main figure uses H_b notation to avoid collisions with metric perturbations. The temperature-only and surface-validation plots remain supplemental.

The baseline article is preserved unchanged. Stage 3 should integrate these results and the Stage 1 bibliography audit, then:

1. Add a compact external-validation paragraph or table with the actual percentage agreements above; do not imply fit-based validation.
2. Present delta M at fixed baryon number = I Omega^2/2 as the general first-law identity for fixed entropy/composition along a regular fixed-EOS sequence. Its roughly 1e-11 monopole residual is a check, not a family-specific discovery.
3. Add Whitman (1985) and the verified primary references; retain the distinction between verified 1982/1983 metadata and papers whose full texts were not obtained.
4. Explain the beta=2 slice of the two-parameter polynomial design and that F(1)=1 is central normalization. A>=0 alone does not certify the whole extended family.
5. Shorten the abstract, round numerical thresholds, state their root-finding origin, and remove repeated limitations.
6. Prepare a physics-journal class variant plus the companion with the user’s requested header, then inspect every compiled page.
7. Prepare code-release metadata and a Zenodo-ready archive. No DOI has been reserved or deposit published; insert a DOI only after a real deposit/reservation.

No new interval certificate was attempted. Existing stability and admissibility thresholds are retained as numerical results. No remaining numerical discrepancy blocks the planned Stage 3 rewrite.

## 5. Reproduction

See the root README for commands. The package includes the original calculation dependencies, benchmark targets, the three extracted source vertices, new calculation scripts and data, seven figure pairs, and the unchanged manuscript baseline. Copyrighted external papers and full external figures are not redistributed. Source hashes and DOI/arXiv provenance are recorded in the data directory.
