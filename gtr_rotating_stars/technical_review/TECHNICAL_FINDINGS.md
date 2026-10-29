# Technical response to the external review

Technical work only. The article source and PDF have not been rewritten.
All reported numerical bounds and thresholds below are floating-point evidence,
not interval-certified theorems. Numerical digits in the data files are retained
for reproduction, not as assertions of that many certified digits.

## 1. Decisions for the next version

The strongest defensible result is a family with an explicit generating
polynomial, exact background quadratures, algebraic first-order dragging and
moment of inertia, and the standard second-order Hartle system with a convergent
integral representation. The calculation now includes tidal deformability,
fixed-baryon mass corrections, fixed-EOS sequences, and radial pulsation estimates.

Three physical limits must be distinguished:

| Property lost | Approximate alpha | M/R | Interpretation |
|---|---:|---:|---|
| Monotone increase of sound speed with density | 0.31508390 | 0.16242453 | Original continuous-stiffening requirement fails above this candidate threshold |
| Positive fundamental radial-mode squared frequency | 0.37628585 | 0.17511372 | Nonrotating configurations become radially unstable near this crossing |
| Causality somewhere inside the star | 8.94592419 | 0.27449835 | Maximum sound speed reaches unity off centre |

These are numerical candidates. The first is obtained from the central sound-speed
slope, supported by radial sampling; it is not a certified global search over the
entire parameter rectangle. The second is an eigenvalue crossing, independently
consistent with a fixed-EOS turning point. The third is a maximization over radius.

**Do not describe the family up to M/R = 0.27 as satisfying all original physical
requirements or as stable stars.** Retain alpha = 0.1 as a transparent benchmark
and alpha = 0.3 as a more compact example satisfying the sampled requirements.
Use the larger-alpha branch as a diagnostic extension with explicit labels.

## 2. What to accept and correct in the review

- Accept the need for family results, a physical interpretation, proper citations,
  conditioning diagnostics, and a conventional introduction and conclusion later.
- The background and alpha = 0.1 rotation values reproduce.
- I–Q agreement is an empirical consistency test, not a verification of Einstein's
  equations and not evidence of novelty.
- The alpha-family changes the EOS. Its mass slope is not a radial-stability test.
- A small matching determinant is not a condition number and depends on basis
  normalization. Scaling the columns resolves much of the apparent problem.
- A rational I for rational polynomial data is a useful consequence of boundary
  matching, not a new general theorem. Arbitrary real polynomial coefficients do
  not imply a rational numerical moment of inertia.
- Volterra/Neumann iteration is standard. Keep the representation and convergence
  argument in an appendix; do not advertise it as a finite closed-form solution.
- The disputed Whitman 1982 reference exists. Reina's 2016 volume and pages are
  also correct. Conversely, the review's Lake and Fodor links are mislabelled.
- Journal suitability and priority cannot be established by the calculations or by
  a short negative literature search.

## 3. Design recipe, with the unexplained coefficients derived

Use R = 1, s = r^2, z = sqrt(1-s), E = 8 pi epsilon, and Pi = 8 pi P.
Let F = b0 + b2 z^2 + b4 z^4 + b6 z^6 + b7 z^7 and prescribe
W = Omega F/N. Define A = 5 F_s + 2s F_ss and D = 2F + s F_s.
The dragging equation gives

    (ln j)_s = -A/(2D),       j(z) = exp[integral_0^z t A(t)/D(t) dt].

For this polynomial,

    A = (-5b2+4b4) + (-14b4+12b6) z^2
        + (35/2)b7 z^3 -27b6 z^4 -35b7 z^5.

Eliminate the constant and z^2 terms, set b7 = alpha, and normalize F(1)=1.
This gives

    b2 = -8alpha/9, b4 = -10alpha/9, b6 = -35alpha/27,
    b0 = 1 + 62alpha/27,
    A = (35alpha/2) z^3(1+2z-2z^2).

The ansatz is a design choice, not a uniqueness claim. It builds the surface
density exponent into the source. F_s is positive on the interval for alpha>0,
so F>=1 and D>0. For example,
35z^4/9-7z^5/2 >= (7/18)z^4.

Write B(r)=integral_0^r F(t)/j(t) dt, with radial composition understood.
The isotropy solution obeys

    d/dr [f/(r^2 j^2 F)] = -2 B/(r^4 j F^2).

At the surface j=1 and j'=0, so

    f'_b = f_b(2+F'_b/F_b)-2B_b/F_b.

Zero pressure requires f'_b=1-f_b. Hence

    f_b = (F_b+2B_b)/(3F_b+F'_b)
        = (F_b+2B_b)/(3N),
    N = F_b+F'_b/3 = 1+26alpha/9.

Matching W to Omega-2J/r^3 then gives the exact family formula

    I/R^3 = 8alpha/[27(1+26alpha/9)].

The surface coefficients are

    E = e3 z^3 + O(z^4),
    Pi = e3 G_b z^5/5 + O(z^6),
    e3 = 35alpha f_b/D(0),       G_b = (M/R)/f_b,
    c_s^2 = G_b z^2/3 + O(z^3).

Consequently P is proportional to n^(5/3) at low density, with a continuous
finite surface. No density threshold or artificial gas cutoff is introduced.

## 4. Family results and I–Love–Q

Q is negative for the oblate configurations. Its tabulated coefficient is
Q/(Omega^2 R^5). The k2 column here is the tidal Love number, not Hartle's metric
function k2(r). Ibar=I/M^3, Qbar=-QM/J^2, Lambda=2 k2/[3(M/R)^5].

| alpha | M/R | I/R^3 | Q/(Omega^2 R^5) | tidal k2 | I–Q fit residual |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 0.08207711 | 0.02298851 | -0.11810929 | 0.12863515 | +0.0792% |
| 0.2 | 0.12896508 | 0.03755869 | -0.09158306 | 0.07774774 | -2.3271% |
| 0.3 | 0.15884974 | 0.04761905 | -0.07786881 | 0.05398075 | -2.8258% |
| 1 | 0.23157988 | 0.07619048 | -0.05600934 | 0.01812094 | -4.9187% |
| 8 | 0.27381404 | 0.09831029 | -0.05311591 | 0.00806212 | -4.9822% |

At alpha=0.1, Ibar=41.57617293, Qbar=18.34360330 and Lambda=23022.79336.
The I–Love fit residual is +1.3448%, not the same as the I–Q residual.
At alpha=0.3 the corresponding I–Love residual is -0.4071%.
The excellent isolated I–Q match does not persist at the 0.1% level along the
family. At alpha=0.03 the polynomial fit gives much larger residuals; this is
low-compactness extrapolation of an empirical fit, not a field-equation failure.

The coefficients used are from Table I of Yagi and Yunes, PRD 88, 023009 (2013):

    ln Ibar = 1.35 + .697 x - .143 x^2 + .0994 x^3 - .0124 x^4,
    x = ln Qbar;
    ln Ibar = 1.47 + .0817 y + .0149 y^2 + .000287 y^3 - .0000364 y^4,
    y = ln Lambda.

The residual definition is Ibar/Ibar_fit - 1. Do not present these approximate
relations as identities of this model.

### Love-number implementation and independent check

Integrate the standard Hinderer variable y=r H'/H with y(0)=2. In s coordinates,

    2s y_s = -y^2 - y F_t - Q_t,
    F_t = [1-s(E-Pi)/2]/f,
    Q_t = {s[5E+9Pi-2E_s/G]/2-6}/f -4s^2 G^2,
    G = Phi'/r = f_s/f-2(ln j)_s.

Use Hinderer's corrected expression for k2 in terms of y(1) and M/R.
There is no surface-density jump because epsilon(R)=0. Independently recover
y(1) from the regular homogeneous static h2,v2 system already used for rotation.
The two values differ by less than 1.9e-11 on the sampled family. This checks
the tidal integration against a different differential system, including the
conversion and surface convention.

## 5. Fixed-baryon mass correction: an exact simplification

Let Mb=m_b N_b, reserving N for the dragging normalization above only in the
code. For a fixed cold EOS the stellar first law is

    dM = Omega dJ + mu_infinity dN_b.

At the nonrotating surface mu_infinity/m_b=sqrt(f_b). Write
M=M0+M2 Omega^2+..., Mb=Mb0+Mb2 Omega^2+..., J=I Omega+....
Along the fixed-central-density rotating sequence,

    M2 = I/2 + sqrt(f_b) Mb2.

Along a fixed-baryon sequence through the same nonrotating star, away from a
turning-point degeneracy,

    delta M|Nb = (I/2) Omega^2 + O(Omega^4).

This is a full-GR consequence of the first law, not a Newtonian kinetic-energy
substitution. The EOS and entropy per baryon must be held fixed. It does not
apply to changes of alpha, which change the EOS.

We independently integrated the baryon correction from the interior:

    rho0 = (epsilon+P) exp(Phi)/sqrt(f_b),
    Mb = 4pi integral_0^R r^2 rho0/sqrt(f) dr,
    delta Mb = 4pi integral_0^R r^2 rho0/sqrt(f)
               [m0/(r f) + r^2 exp(-2Phi) W^2/3 + p0*/c_s^2] dr.

All three bracketed terms have order Omega^2. Surface-motion terms vanish
because rho0 vanishes at the surface. The apparently divergent p0*/c_s^2 is
integrable after multiplication by rho0.

| alpha | M2/R^3, fixed centre | Mb2/R^3, fixed centre | M2/R^3, fixed baryons |
|---:|---:|---:|---:|
| 0.1 | 0.3562572253 | 0.3771008761 | 0.01149425287 |
| 0.3 | 0.2514770015 | 0.2756215333 | 0.02380952381 |

The first-law residual is below 2.4e-11 in the dimensionless coefficient across
the nine sampled models. This is a stronger new internal check than the empirical
I–Q fit. To construct all fixed-baryon monopole functions, add the fixed-EOS
background derivative with central-density shift

    delta E_c/Omega^2 = (I/2-M2)/(dM/dE_c)_EOS.

It is approximately -8.48205 at alpha=0.1 and -88.1636 at alpha=0.3 in R=1
units. This is not an alpha derivative. The shift grows near a radial turning
point, where the fixed-baryon expansion needs separate treatment. The quadrupole
and I at the orders retained are unchanged by this O(Omega^2) background shift.

## 6. Radial stability: two independent calculations

For each reference alpha, reconstruct one EOS on its occupied density range,
then integrate new TOV configurations with lower central enthalpy using that
same EOS. The table uses h_c/h_c,reference from 0.002 to 1; it never invokes the
arbitrary high-density completion. The EOS is interpolated using E^(2/3) and
Pi^(2/5) as functions of enthalpy to regularize their gas-surface powers.

Hydrostatic pressure is reconstructed by an integral near the surface. Direct
subtraction of terms of order z^3 to obtain pressure of order z^5 otherwise
produces spurious minima of Gamma. This numerical problem was found and fixed.

The reference star is recovered with |Delta R| and |Delta M| below about 6e-11
for alpha<=1, and about 2.5e-9 and 3.6e-10 respectively at alpha=8.

We also solve the nonrotating radial Sturm–Liouville problem by Rayleigh–Ritz,
using zeta=r^2 exp(-Phi) xi and trial functions r^3 P_k(2r^2-1):

    (Pi_r zeta')' + (Q_r + sigma^2 W_r) zeta = 0,
    Pi_r = exp(Lambda+3Phi) Gamma P/r^2,
    Q_r = exp(Lambda+3Phi)/r^2
          [(P')^2/(epsilon+P) -4P'/r -8pi P(epsilon+P) exp(2Lambda)],
    W_r = (epsilon+P) exp(3Lambda+Phi)/r^2,
    exp(2Lambda)=1/f, Gamma=(epsilon+P)c_s^2/P.

Here Pi_r,Q_r,W_r are pulsation coefficients, not the EOS pressure or stellar
quadrupole. The centre condition is zeta=O(r^3), and the natural zero-pressure
surface condition is Pi_r zeta' ->0. The equilibrium EOS derivative is used as
the adiabatic derivative, consistently with the assumed isentropic fluid.

| reference alpha | dM/dE_c at fixed EOS | sigma0^2 R^2 | Interpretation |
|---:|---:|---:|---|
| 0.1 | +0.0406462 | +0.24392502 | Numerically stable |
| 0.2 | +0.0110838 | +0.17425835 | Numerically stable |
| 0.3 | +0.00258233 | +0.07153945 | Numerically stable |
| 0.31508390 | +0.00191039 | +0.05669432 | Numerically stable |
| 1 | -0.00199874 | -0.31746788 | Radially unstable |
| 8 | -0.00006653 | approximately -0.576 | Radially unstable |

E_c=8pi epsilon_c in this table. One-sided derivative step convergence is saved
separately. With 8, 16, and 24 trial functions, the low-alpha fundamental
eigenvalues agree to better than 2e-9 between the smallest and largest bases;
the 16-to-24 change is a few 1e-12. Higher-alpha modes converge more slowly;
only their instability sign is needed here.

The zero crossing is alpha≈0.37628585. At that alpha, a separate fixed-EOS
mass derivative approaches zero as the finite-difference step shrinks.
Positive Ritz values are upper bounds, not rigorous positive lower bounds for
the true eigenvalue. These computations provide converged numerical stability
evidence, supported by the independent sequence calculation; they are not a
validated spectral proof. Rotational corrections to the radial frequencies and
nonaxisymmetric stability have not been calculated.

## 7. Energy conditions, scales, and conditioning

The sampled backgrounds have positive density and pressure, decreasing density,
and P/epsilon<1. This implies the weak, dominant and strong energy conditions
for the perfect fluid on those samples. At alpha=0.1, max(P/epsilon)≈0.05286;
at alpha=0.3 it is about 0.1281. The causal-threshold model has max(P/epsilon)
about 0.851, so loss of DEC is not what sets the reported causal threshold.
An interval certificate is still needed to replace sampling with a global proof.

The Buchdahl condition is M/R<4/9 under its usual static isotropic monotone-density
assumptions. Even the broad causal branch reaches only approximately 0.2745;
the original monotone-stiffening branch reaches approximately 0.1624.

At R=12 km the alpha=0.1 and 0.3 models have masses approximately 0.6670 and
1.2909 solar masses. The characteristic frequency

    f_dyn = (c/2pi R) sqrt(M/R)

is approximately 1139 and 1585 Hz respectively. It equals the Schwarzschild
circular-geodesic frequency at the unperturbed surface, as measured at infinity.
It is a dynamical scale, not a computed rotating mass-shedding limit. Require
eta=Omega^2 R^3/M <<1, and separately check the actual metric and surface
perturbations, especially when holding baryons fixed near a turning point.
For example eta=.01 corresponds to about 114 and 158 Hz. No accuracy percentage
or actual Kepler limit follows without higher-order or rapid-rotation calculations.

At alpha=.1 the unscaled quadrupole matching matrix has 2-norm condition number
about 2661, reduced to 12.47 after unit-norm column scaling. Its small determinant
mostly reflects the normalization of the decaying exterior solution. The code
now scales these columns before solving and evaluates Q_2^2 and the related
vacuum basis by convergent inverse-argument series, avoiding logarithmic
cancellation at small compactness. This basis representation converges for
t=r/M-1>1; the 80 terms used are ample for the sampled t>=2.6.

## 8. Tightening the existing small-alpha admissibility argument

The earlier continuity argument needs endpoint justification, not merely a claim
of uniform continuity. It can be completed locally in parameter as follows.

1. D>0 gives analytic parameter dependence of j's rational integral on the
   compact physical interval for alpha in a neighbourhood of zero.
2. At the centre use u=(f-1)/s. The regular Frobenius solution is analytic in s
   and alpha. The boundary equation for its free central coefficient has
   denominator 3h(1)+2h_s(1)=3 at alpha=0, so this coefficient is analytic in a
   sufficiently small parameter neighbourhood.
3. At the surface use z. The quadratures have analytic coefficients there;
   A is exactly divisible by z^3 and the enforced Pi(0)=0 gives Pi=z^5 times a
   regular function. E=z^3 times a regular function follows from E+Pi=-4fa.
4. Match these endpoint charts to the nonsingular interior. The functions
   E/(alpha z^3), E_z/(alpha z^2), Pi/(alpha^2 z^5),
   c_s^2/(alpha z^2), and (c_s^2)_z/(alpha z), with their endpoint and alpha=0
   limits, are jointly continuous on a compact rectangle for sufficiently small
   alpha. Their leading coefficients are strictly positive, as in the original
   Bernstein-polynomial argument. Compactness now gives the uniform sign result.

This proves existence of some admissible nonempty interval; it still does not
give a numerical upper bound or certify alpha=.1. It must not be replaced by
the unproved assertion that all sampled values form a certified interval.

An explicit certificate should target a conservative closed interval, e.g.
[0.05,0.30], then handle the small-alpha interval separately. Enclose the regular
centre coefficient and the background using validated quadrature or a validated
ODE solver, propagate enclosures through the factored expressions above, and
check f>0 and 1-c_s^2>0 as well. The endpoint Taylor remainder and derivative
enclosures are essential. Merely calling Arb on sampled quadratures or enclosing
floating-point output does not certify the inequalities. No such machine
certificate was completed in this pass; python-flint is not installed in the
current runtime. No explicit certified alpha interval is asserted.

The C-infinity high-density completion previously constructed applies when the
central stiffening derivative is strictly positive. Use a strict interior
parameter range, not the exact monotonicity threshold, for that argument.

## 9. Remaining technical work, ordered by importance

| Priority | Task | Acceptance condition |
|---|---|---|
| P0 | Resolve the relation to Whitman 1982/1983 | Inspect full texts and compare the exact differential identity and inverse choice of dragging profile; distinguish classical method from specific new ansatz |
| P0 if a rigorous admissibility claim is retained | Validated interval certificate | Endpoint-factored global inequalities with explicit error enclosures, preferably on [0.05,0.30]; otherwise say numerical admissibility |
| P1 | Matched-compactness comparison | Recompute Tolman VII, a specified Buchdahl model and a thermodynamically consistent n=1 polytrope with the same I/Q/Love conventions; do not import incomparable table values |
| P1 | Final figure suite | Add EOS comparison and a surface-shape figure with an explicit gauge/geometric-radius definition; diagnostic profile and family figures already exist |
| P1 | Literature-comparison completion | Inspect Stephani chapter 16 and relevant rotating-fluid sections, and compare the specific invariants/EOS to catalogued solutions |
| Optional stronger validation | Validated radial spectrum or separate shooting solver | Positive lower bounds for sigma0^2 if claiming rigorous stability; existing converged Ritz results suffice for a clearly numerical claim |
| Optional extended physics | Rotation-induced frequency shifts and actual mass shedding | Separate calculation; not inferred from the second-order metric or f_dyn alone |

The finite standard-special-function solution of the second-order interior
remains unavailable. There is no new proof of nonintegrability. Neither a
literature audit nor additional numerical physics removes that obstruction.

## 10. Rewrite plan for the next prompt only

1. Preserve the supplied IEEEtran/fallback header, changing the title only as
   appropriate. Establish one notation table first: Newton's G, g=Phi'/r,
   pressure P, Legendre polynomial P_2, tidal k2 versus metric k2, baryon number,
   normalization constant, and background/perturbation coefficients.
2. Introduce the exact-family benchmark objective and existing generating
   literature; avoid novelty claims not cleared by the full-text comparison.
3. Show the polynomial recipe and boundary normalization derivation before
   presenting the quadratures. Explain the finite adiabatic gas surface once.
4. Present the parameter family with the three separate physical thresholds.
   Lead physical examples with alpha=.1 and .3, not the unstable causal endpoint.
5. Present Hartle matching, I, Q, Love number, baryon convention, radial
   stability, and physical frequency scales with common conventions and figures.
6. Move the Volterra representation, endpoint proof, continuation details, and
   reproduction instructions to appendices/supplement. State the remaining
   finite-closed-form limitation once.
7. Conclude with demonstrated results and their scope. Code filenames belong
   in a supplement, not in the scientific narrative. Journal selection comes
   after the novelty comparison and remaining technical decisions.

## 11. Reproduction and files

Run the supplied Python files with numpy, scipy and matplotlib installed:

    python technical_family.py
    python technical_stability.py
    python technical_radial.py
    python technical_verification.py
    python technical_figures.py

These produce family tables, threshold candidates, fixed-EOS sequences,
radial-mode convergence data, Einstein residual checks and four diagnostic
figures in technical_review/. `drag_regular.py`, `drag_second.py` and
`check_star.py` are included dependencies. No external stellar-EOS tables are
used. See REFERENCES_AND_PRIORITY.md for the source audit.

The Einstein residual check is a finite-difference diagnostic at representative
interior points. The diagonal residual scales as Omega^4 and the t-phi residual
as Omega^3 above its differentiation noise floor, as expected for the retained
metric. This supports, but does not replace, the analytical reduction to the
Hartle equations. The exact background quadratures and algebraic first-order
dragging retain full GR dependence throughout; no weak-field background is used.
