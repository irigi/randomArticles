# Rotating-star revision: Stage 1 audit and handoff

16 October 2026. Stage 1 of the requested three-stage revision is complete.
The article itself has not yet been rewritten in this pass. Stage 2 is the
external numerical validation and figure preparation; Stage 3 is the integrated
article, bibliography, PDF and submission/source package.

## 1. Supplied files and baseline

Both supplied ZIPs were extracted into separate directories and inspected.
The corrected Figure 1 archive retains the original background, rotation,
stability and comparison calculations and data. Among corresponding calculation
files, only the figure-generation script differs. The thermodynamics archive
uses a byte-identical copy of the background solver.

The corrected Figure 1 resolves the thin outer layer by sampling closer to the
surface, computes small pressures from the hydrostatic integral, and marks the
analytic endpoint. Its shortened density/pressure axis label is visible.
These changes will be carried into the revision. The supplied article also
clarifies that the 5/3 law is a surface asymptote, not a global polytrope; this
wording will be retained. Its document-class change will be handled separately
from the physics changes.

## 2. Whitman papers: what is actually established

The supplied six-page paper is P. G. Whitman, *Interior solutions for rotating
fluid spheres*, Physical Review D **32**, 1857–1862 (1985),
doi:[10.1103/PhysRevD.32.1857](https://doi.org/10.1103/PhysRevD.32.1857).
Its text and equations were examined, including the first two rendered pages.

It studies analytical first-order dragging on static spherical backgrounds.
Sections II–III transform the dragging equation into Heun or hypergeometric
form; Section IV develops matched dragging and moment-of-inertia examples.
It does not supply the general second-order monopole and quadrupole construction
required in the present paper. In its notation, its metric function f is our
j², its u is our s, and its relative angular velocity has the opposite sign
to W. Those correspondences make the connection direct.

Its bibliography explicitly contains both disputed citations:

- Reference 11: Whitman, Physics Letters **89A**, 129 (1982).
- Reference 17: Whitman, Physical Review D **27**, 1722 (1983).

Publisher-deposited DOI metadata independently confirms their titles and page
ranges: *The Tolman-Wyman equation in general relativity*, Physics Letters A
**89**, 129–130 (1982), doi:[10.1016/0375-9601(82)90873-8](https://doi.org/10.1016/0375-9601(82)90873-8);
and *Solutions to the general-relativistic Tolman-Wyman equation*, Physical
Review D **27**, 1722–1727 (1983),
doi:[10.1103/PhysRevD.27.1722](https://doi.org/10.1103/PhysRevD.27.1722).
The APS abstract for the latter is accessible and describes its static
generalization procedure.

**Access limitation:** the complete 1982 and 1983 papers were not obtained in
this pass. Their existence and bibliographic identities are confirmed; their
full derivations must not be described as independently inspected. The 1985
paper supplies a directly inspected primary source for the analytical-dragging
context. The revised introduction will cite it and will make no priority claim
for the generating strategy or for the specific polynomial family.

## 3. External quadrupole benchmark selected

Two complementary routes are now available. Agreement has not yet been tested;
that is the central Stage 2 calculation.

### Published numerical table

Chandrasekhar and Miller, *On slowly rotating homogeneous masses in general
relativity*, MNRAS **167**, 63–80 (1974),
doi:[10.1093/mnras/167.1.63](https://doi.org/10.1093/mnras/167.1.63),
Table I on p.73, gives explicit I and Q values. The scanned table and its
units were visually inspected. With R_S=2M, its columns are I/R_S³ and the
positive oblate quadrupole in units J²/R_S. Thus our Ibar is eight times its
I column and our Qbar is half its Q column.

| R/R_S | C=M/R | Published I/R_S³ | Published Q column | Derived Ibar | Derived Qbar |
|---:|---:|---:|---:|---:|---:|
| 2 | 0.25 | 1.049 | 6.525 | 8.392 | 3.2625 |
| 3 | 1/6 | 2.123 | 15.67 | 16.984 | 7.835 |
| 5 | 0.10 | 5.487 | 37.73 | 43.896 | 18.865 |

These are rounded historical numerical values, not exact targets. Stage 2
will test the comparison solver using the analytic homogeneous background,
with the correct nonzero-surface-density treatment. This is a code benchmark;
it does not make an incompressible EOS a candidate physical model for our star.

### Published polytropic curves

The original arXiv source archive for Yagi and Yunes, PRD **88**, 023009 (2013),
was obtained from [arXiv:1303.1528](https://arxiv.org/abs/1303.1528).
It includes the vector files IbarC-PRD.eps and QbarC-PRD.eps. The n=1 curves
can therefore be digitized with a documented coordinate precision, rather
than compared with the universal fit. This is a separate published result.

**EOS convention matters:** their equations (6) and (9) use p=K rho² with
rho the total energy density. Our existing comparison polytrope instead uses
P=K rho0², epsilon=rho0+P. They are different relativistic EOSs. Stage 2 must
implement the former for this particular benchmark and retain the latter in
the current equal-compactness table. A figure-derived comparison must state
its digitization uncertainty; it cannot claim the precision of the ODE solver.

Urbanec, Miller and Stuchlík (2013),
[arXiv:1301.5925](https://arxiv.org/abs/1301.5925), primarily supplies realistic
and strange-star EOS curves, so it is less convenient for this exact polytropic
test. The existence of several attractive references is not itself a completed
validation.

## 4. Reference-model I–Q residuals

Recomputing the residuals directly from the archived full-precision comparison
results gives the following percentages, using Ibar/Ibar_fit(Qbar)-1:

| Model | C=0.0820771141 | C=0.1588497390 |
|---|---:|---:|
| Present family | +0.0792% | -2.8258% |
| Tolman VII | +1.5892% | -0.5387% |
| Rest-mass index-one polytrope | -0.1410% | -2.3089% |
| Buchdahl gas | +1.2130% | -0.1033% |

Stage 2 will add these reference-model points to the I–Q residual panel and
calculate their I–Love residuals consistently. The original family curve and
the distinction between empirical-fit comparison and correctness testing will
remain clear.

## 5. Thermodynamic integration and necessary wording corrections

The supplement is suitable for a new four-panel figure after the EOS
reconstruction discussion. Its quantities are temperature-equivalent P/n,
enthalpy per particle, internal energy above rest mass, and the effective
polytropic coefficient. The separate temperature-only plot is redundant in
the main article; it can remain in the supplement. The surface-Gamma zoom is
useful as a supplementary validation figure.

The following analytical checks were completed independently from the surface
expansion. With z=sqrt(1-r²/R²), g_b=C/(1-2C), and H=epsilon+P,

    H = H3 z³ [1+2z+O(z²)],
    P = (H3 g_b/5) z⁵ [1+(5/3)z+O(z²)],
    n = n3 z³ [1+2z+O(z²)],
    Gamma = 5/3 - (5/9)z + O(z²).

Thus the final approach to 5/3 is from below. The shallow dip should be
explained, not removed or smoothed away. Stage 2 will reproduce the supplied
near-surface numerical checks and quantify the plotting error.

Write the enthalpy per particle as H_b=(epsilon+P)/n, with a distinct typeset
symbol in the article. Then

    H_b/m_b = sqrt(f_b) exp(-Phi),
    theta = P/(m_b n) = (g_b/5) z² [1-z/3+O(z²)],
    u_int = epsilon/(m_b n)-1 = H_b/m_b-1-theta,
    u_int/theta -> 3/2.

In units c=1, theta equals k_B T_ideal/m_b only under the extra identification
P=n k_B T. The isentropic barotrope fixes neither a unique physical temperature
nor heat capacity: those require epsilon(n, entropy per particle). This
limitation will be stated once beside the new figure. The caption should use
“ideal-gas-equivalent temperature,” rather than imply a fully specified thermal
EOS throughout the core.

Two existing sentences need correction:

1. H_b tends to m_b at the surface. It is the rest-mass-subtracted enthalpy,
   and the logarithmic enthalpy, that vanish.
2. H_b is not generally the Gibbs chemical potential at finite entropy.
   The notation should say enthalpy per particle; the chemical potential is
   H_b-T times entropy per particle. This also clarifies the first-law section.

The star is in mechanical equilibrium. A Tolman condition for complete thermal
equilibrium is not being imposed on the conditional temperature profile.

## 6. Fixed-baryon mass correction

Accept the review's correction in substance. At fixed EOS, composition and
entropy per particle, along a regular uniformly rotating equilibrium sequence
with fixed baryon number, the general first law gives

    dM = Omega dJ,
    J = I0 Omega + O(Omega³),
    M(Omega)-M(0) = (1/2) I0 Omega² + O(Omega⁴).

This is a general identity, not a special property of this family. Its agreement
with the independently computed monopole and baryon corrections is the useful
result here. The qualifier concerning fixed entropy/composition should remain;
“every star” without those qualifications is too broad.

Bardeen (1970), Section IX, equations (79)–(80), supplies the primary source.
Equation (80) explicitly handles variations at fixed entropy per baryon and
discusses their use in checking numerical calculations. Replace the lecture
notes citation with this paper. The redshifted coefficient multiplying dN on
an isentropic sequence is enthalpy per particle, not an unqualified finite-
temperature chemical potential.

## 7. The free polynomial coefficient

The review is correct that the generating ansatz has another free parameter.
One convenient parametrization is b7=alpha and b6=-35 beta alpha/54. Eliminating
the constant and z² terms in A and imposing central normalization gives

    F(z)=1-alpha+(89 beta alpha/54)
         -(4 beta alpha/9)z²-(5 beta alpha/9)z⁴
         -(35 beta alpha/54)z⁶+alpha z⁷,
    A=(35 alpha/2)z³(1+beta z-2z²).

For alpha>0, beta>=1 makes A>=0 on [0,1]. This is a necessary design property,
not a proof that every such background is physically admissible. At beta=1,
A vanishes at the centre and epsilon_c+P_c=0; that endpoint cannot be a
positive-density, nonnegative-pressure central fluid state. The published slice
is beta=2. Its positive factorization, known admissibility and explicit inertia
make it a sufficient focused example; a complete two-parameter survey is not
needed for this revision.

The generalized normalization and inertia are

    N=1-alpha+35 beta alpha/18,
    I/R³=4 beta alpha/(27 N).

The article will state that F_polynomial(1)=1 is normalization at r=0 because
z=1 there. N, obtained by vacuum matching, fixes the physical normalization at
infinity. This avoids confusing central and surface values.

## 8. Bibliography changes prepared

The following primary references are identified, with DOI metadata or full
papers checked as indicated:

- Bardeen, *A Variational Principle for Rotating Stars in General Relativity*,
  ApJ **162**, 71–95 (1970), doi:10.1086/150635. Full paper inspected for the
  first-law and entropy conventions, especially equations (79)–(80).
- Chandrasekhar, *The Dynamical Instability of Gaseous Masses Approaching the
  Schwarzschild Limit in General Relativity*, ApJ **140**, 417–433 (1964),
  doi:10.1086/147938. Original paper obtained. Use for radial pulsations,
  together with Bardeen–Thorne–Meltzer's explicit Sturm–Liouville formulation.
- Bardeen, Thorne and Meltzer, *A Catalogue of Methods for Studying the Normal
  Modes of Radial Pulsation of General-Relativistic Stellar Models*, ApJ
  **145**, 505–513 (1966), doi:10.1086/148791. Primary text accessed.
- Hinderer, erratum to *Tidal Love Numbers of Neutron Stars*, ApJ **697**,
  964 (2009), doi:10.1088/0004-637X/697/1/964. Publisher-deposited bibliographic
  metadata checked. The corrected arXiv v4 explicitly identifies corrections
  to equations (20), (23) and Table 1; the separate publisher erratum PDF
  could not be fetched. Cite the original paper and the published erratum.
- Yagi and Yunes, *I-Love-Q: Unexpected Universal Relations for Neutron Stars
  and Quark Stars*, Science **341**, 365–368 (2013),
  doi:10.1126/science.1236462, arXiv:1302.4499. Add alongside the detailed PRD
  paper, which remains the source of the fit coefficients and polytropic curves.
- Jiang and Yagi, *Analytic I-Love-C relations for realistic neutron stars*,
  PRD **101**, 124006 (2020), doi:10.1103/PhysRevD.101.124006,
  arXiv:2003.10498. Their compactness-expanded Tolman VII work belongs in the
  comparison context; it is not an exact finite-compactness replacement for
  our background quadratures.
- Whitman (1985), as discussed above, is added to the historical context.

## 9. Remaining editorial and delivery decisions

- Shorten the abstract and round threshold compactnesses to four significant
  figures. Describe them as numerical boundaries.
- Explain that thresholds were obtained by root-finding and radial maximization,
  not read off the few points in the family plot.
- Trim the final caveat paragraph; keep each material limitation where needed.
- Describe the surface embedding as an intrinsic geometry of a specified
  constant-time section, not a foliation-independent stellar shape observable.
- Prepare a physics-journal RevTeX version in Stage 3, retaining a companion
  version with the user's requested header. RevTeX was not found in the current
  TeX installation; compilation support will need to be arranged then.
- Prepare a Zenodo-ready code/data release, metadata and citation file in
  Stage 3. A DOI can be inserted only after an actual deposit or reservation.
  No publication, account action or DOI reservation has occurred in Stage 1.
- Interval certification is optional follow-up, not a condition for completing
  this three-stage revision. Numerical boundaries must remain labelled as such.

## 10. Stage 2 execution order

1. Implement the homogeneous-background response test and compare the published
   I and Q table values with explicit sign/normalization conversions.
2. Decode the published Yagi–Yunes vector curves and their coordinate precision;
   implement their energy-density n=1 EOS and compare at identical compactness.
3. If either test disagrees beyond source precision, investigate the discrepancy
   before revising the article's validation claim.
4. Add the existing reference-model I–Q/I–Love residuals to the universal-relation
   figure without treating those fits as an independent validation.
5. Reproduce the thermodynamic diagnostics, surface Gamma law and corrected
   Figure 1 with consistent notation and visual style. Record results for Stage 3.

Stage 3 then integrates only the results actually established, compiles the
article, checks every page, and packages the sources, figures and reproduction
instructions. This checkpoint intentionally stops before that next calculation
stage, as requested.
