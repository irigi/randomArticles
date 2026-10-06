# Plan — Osmosis lab: chemical potential with mobile ring hosts (draft, 2026-10-06)

Status: proposal, not started. Supersedes the "hosts inside the Carnot cylinder"
draft. Source: the original brief, §8.3 and §10 (esp. §10.6).

## 0. Carnot is kept as it is, behind its preset

Nothing changes for Carnot: `--preset carnot_discs` and `--preset carnot_triangles`
already select it in `run`, `precalculate`, `gui` and `replay`, and the Carnot-only
options are already checked against the `carnot_` prefix. The new
demonstration uses its own preset family:

    microthermo gui --preset osmosis            # ring hosts with a binding well
    microthermo gui --preset osmosis_hard       # same rings, ε = 0 (control)

with `osmosis_`-only options (`--hosts`, `--binding-energy` (default 1.5),
`--host-mass`, `--host-area-fraction`, `--mouth-width`, `--rough-fraction`),
validated the same way. The three old stub presets `labyrinth`,
`labyrinth_energetic` and `selective_membrane` are removed, because their fixed walls and fixed portals are the model we are
moving away from. The window picks `LabView` for `carnot_*` and a new
`OsmosisView` for `osmosis*`. Carnot stays the default until the osmosis gates
in §8 pass.

## 1. Why this demonstration

Chemical potential is defined operationally by **what stays equal when particles
can move between two systems and nothing else can**. A semipermeable membrane
makes that literal:

- The **right chamber** holds only small discs. It is a plain ideal gas, so its
  density gives μ directly: μ/T = ln(c_R) + const. It works as a **μ reference
  electrode**.
- The **left chamber** holds the same discs plus mobile ring hosts that can
  capture them. The membrane lets discs through and stops hosts.
- At equilibrium, μ_D is equal on both sides. Concentration, total disc count
  and pressure are not.

Every headline below comes from trajectories and wall impulses. None is put
into the dynamics.

1. **Equal μ, unequal counts.** The free-disc density on the left (away from
   hosts) equals c_R. The total disc count on the left differs, and the sign
   depends on ε/T.
2. **The crossover with temperature (headline).** A host removes area from the
   free gas (its footprint) and returns a cavity weighted by e^{ε/T}. For the
   default geometry the crossover is near ε/T ≈ 1.2. With ε = 1.5, **cooling
   pulls discs through the membrane into the host side, and heating pushes
   them back out**, past parity. With ε = 0 nothing happens on heating (the
   hard-only control).
3. **Osmotic pressure counts the hosts, not the discs.** The net force on the
   membrane is Π ≈ n_hosts·T (van 't Hoff), however many discs are trapped,
   because the discs are equally free on both sides. Measured from impulses on
   the membrane.
4. **Mouth width changes the rate, not the equilibrium** (kinetics). Narrow and
   wide mouths with the same cavity and ε reach the same equilibrium; only the
   relaxation time differs.

## 2. Bodies and geometry

All bodies are circles or pieces of circles. Every time-of-impact (TOI) is
then a quadratic with an exact solution (§5).

**Small gas discs D:** mass m = 1, radius r (default 0.02), with spin ω_d and
moment of inertia I_d = ½·m·r² (uniform disc). The spin only changes in rough
contacts with hosts (§3).

**Ring hosts L:** one rigid body each, **mobile and rotating**: position,
velocity, angle θ, angular velocity ω, mass M (default 25), moment of inertia
I = ½·M·(R_in² + R_out²) (uniform annulus). The mass is **declared balanced**:
the centre of mass sits at the geometric centre, as if a small counterweight
opposite the mouth made up for the missing wall. This is an explicit
configuration, recorded in the metadata, as the brief requires for incomplete
rings. The geometry in host coordinates:

- an annular wall from R_in to R_out, with wall thickness h = R_out − R_in;
- a **mouth**: a gap of clear chord width w (default 3 disc diameters) at body
  angle 0; the initial angles θ are random;
- two **end caps** closing the wall at the mouth: circles of radius h/2 centred
  at the mid-wall radius;
- a **binding well**: a step circle of radius R_s = R_in − 2r, concentric with the
  host. A disc centre inside it has U = −ε. The thin neutral band between R_s
  and the wall avoids simultaneous wall-and-step events.

Collision features: outer arc (contact only when the contact angle lies outside
the mouth), inner arc (likewise), the two end caps, and the step circle (not
solid; see §3).

### 2.1 How rotation exchanges energy, cheaply

**The problem.** A frictionless circular arc cannot be spun: every contact
normal passes through the host centre, so the impulse has no torque. Only the
two small end caps would give torque, so the rotation would hardly couple to
the gas.

**The choice: rough contacts on circular hosts** (recommended). Keep the ring
geometry and change the contact law. Every host contact (disc–host,
host–host, host–post) is **perfectly rough**: the *whole* relative velocity of
the two contact points is reversed, normal and tangential components alike.
This is the classical rough-sphere model of kinetic theory (Bryan 1894,
Pidduck 1922; Chapman & Cowling, ch. 11). It is:

- **Energy-conserving and elastic.** Reversing the contact velocity is a
  reflection in the mass metric, so kinetic energy is unchanged.
- **Momentum- and angular-momentum-conserving.** The impulses are equal and
  opposite at the same contact point. This is why the discs carry spin: a
  spinless point disc would break angular momentum by r × j_t.
- **Time-reversible and measure-preserving.** The Maxwell–Boltzmann
  distribution, including the rotational modes, stays invariant, so the
  equilibrium statistics and the μ references are unchanged.
- **Strongly coupling.** A tangential impulse at the contact point spins the
  host on almost every contact, so host rotation thermalises on the time scale
  of collisions.

Formula: with G the 2×dof contact Jacobian (normal and tangential rows) and
D = G·M⁻¹·Gᵀ (2×2), the impulse is j = −2·D⁻¹·g and the velocities change by
u⁺ = u⁻ + M⁻¹·Gᵀ·j. On an arc, the normal passes through the host centre, so
D is diagonal. The rule is then two independent 1-D reflections (normal and
tangential), each with a scalar formula. On an end cap, D has an off-diagonal
term; the 2×2 inverse is still trivial.

**Tunable coupling (optional).** With probability p_rough the contact is rough,
otherwise smooth. Both maps preserve energy and the equilibrium distribution, so any
mixture does too. The random draw comes from the same uniform buffer as the
thermal walls. Default p_rough = 1; lowering it slows rotational relaxation
without changing equilibrium, a nice "kinetics, not equilibrium" knob.

Disc–disc contacts stay smooth (cheaper, and the spin bath reaches the gas
through the hosts). Thermal walls refresh only the normal mode, as now.

**Why this is cheap: rotation does not move a circle.** The outer arc, inner
arc and step circle are all centred on the host centre, so their positions do
not depend on θ. Their TOI is the same quadratic as without rotation. Rotation
enters only through two things:

1. **Validity check at the root.** A root on an arc counts only if the contact
   angle lies outside the mouth at that time, θ(t*) = θ₀ + ω·t*. This is one
   evaluation, done once per root.
2. **End caps, inside a bounded "mouth window".** A cap's centre always lies on
   the mid-wall circle, so the cap never leaves the wall band
   [R_in, R_out]. A disc can therefore touch a cap only while its centre is in
   the band [R_in − r, R_out + r]. The entry and exit times of that band are a
   quadratic (radial distance does not depend on θ). It can enter the band
   only through the mouth, since elsewhere it hits an arc first. Only inside
   that short interval is a cap TOI searched, by conservative advancement:
   the distance changes at most at |v_rel| + |ω|·R_mid, then a bracketed root.
   This uses the same pattern as `ccd.py`. Discs cross the mouth rarely, so
   this costs little. The same window applies to host–host, host–post and
   host–wall contacts: caps never reach past R_out, so they can touch something
   only where that something is inside a mouth.

**Alternatives considered:**

| Option | How rotation couples | TOI cost | Verdict |
|---|---|---|---|
| **Rough circular ring** | tangential impulse on every contact | quadratic + rare mouth-window search | **recommended** |
| Smooth ring, centre of mass off the geometric centre | normal impulses gain a lever arm e | the circle centre wobbles with θ, so *every* host TOI becomes an iterative search; coupling ∝ e is weak | no |
| Smooth polygon with a gap (triangle, pentagon) | corners and faces give torque naturally | every face rotates, so every host contact needs conservative advancement on non-convex geometry (the slow triangle path, but worse) | no; most expensive |
| Smooth ring, I = ∞ (earlier draft) | caps only, no rotation | all quadratics | superseded by your request |

The polygon is the "purest" frictionless model, but it pays the price of the
current triangle engine on every host contact. The rough ring keeps exact
quadratics for almost all events and needs no non-convex distance queries.

**Why mobile, rotating hosts matter.** A disc trapped in a fixed circular
cavity is an integrable billiard: its angular momentum about the centre is
conserved, so its radial speed at the step never changes. It escapes now or
never. With rough contacts, every hit on the inner wall exchanges angular
momentum and energy with the spinning host, and kicks from outside discs move
the wall under the trapped disc. Trapped discs therefore gain or lose the
energy needed to climb the step, and the step impulse feeds binding energy into
host motion. Gate E0a checks this mixing: with a frozen host the disc must
never escape, and with a free host it must.

**One ring, not concentric rings, at first.** A single ring with one mouth
gives all four headlines. A second inner ring with its mouth turned away is
built from the same primitives (more arcs on the same body), so it is a cheap
later variant for headline 4 if the mouth-width contrast is too weak.

**Membrane:** a column of fixed circular posts at x = x_m.

- Gap between posts > 2r + clearance, so discs pass.
- Post diameter > w, so a post can never enter a host mouth.
- Gap < 2R_out, so a host can never pass.

These three are checked as geometry tests at build time. The membrane is
insulating (specular).

**Box:** a rectangle of 4 × 2, membrane at x = 2. Outer walls are thermal at T
for all bodies (the existing flux-weighted kernel works for any mass), so
hosts are thermalised directly as well as through discs.

**Host sizing as their number grows (your point 2).** The hosts' footprint
fraction φ of the left chamber is held fixed (default 0.2), so
R_out = √(φ·A_L/(π·N_L)). Wall thickness and mouth width are tied to the
**disc** size, not the host size, because they must stay passable and
impenetrable at any scale. The cavity therefore shrinks faster than the host:

| hosts | R_out | R_in | step R_s | cavity a | crossover ε*/T |
|---:|---:|---:|---:|---:|---:|
| 4 | 0.252 | 0.222 | 0.182 | 0.104 | 0.80 |
| **8** | **0.178** | **0.148** | **0.108** | **0.037** | **1.21** |
| 16 | 0.126 | 0.096 | 0.056 | 0.010 | 1.91 |
| 32 | 0.089 | 0.059 | 0.019 | 0.001 | 3.48 |

(r = 0.02, h = 0.03, φ = 0.2.) Above about 16 hosts the cavity holds less than one
disc. Going further needs smaller discs, scaled like the Carnot
`radius ∝ 1/√N` rule. The builder enforces R_s ≥ 3·(2r) and reports the
limit, rather than building a useless host.

## 3. Event physics (all reused or derived from the brief)

- **Disc–disc, disc–wall, disc–post:** smooth hard elastic impulses, as now.
- **Host–disc, host–host, host–post:** perfectly rough (§2.1),
  j = −2·D⁻¹·g, with probability p_rough (default 1), otherwise smooth.
  Host–wall contacts use the wall's own policy (specular or thermal normal
  refresh); host spin is unchanged there.
- **Step crossing:** this is the brief's §8.3 law with a moving host.
  g = relative radial speed, D = 1/m + 1/M, and the normal passes through the
  host centre, so there is no torque even with rotation.
  - A crossing is allowed if g² ≥ 2·D·ΔU. Then g' = sign(g)·√(g² − 2·D·ΔU) and
    j = (g' − g)/D, applied equal and opposite to the disc and the host.
  - Otherwise the disc reflects: g' = −g.
  - The membership label flips only on a crossing.
  - This generalises the existing `resolve_energy_step`; its threshold
    convention is kept.
- **Thermal walls:** the existing flux-weighted normal-mode refresh. Q is
  logged per wall and per species.
- **Ledger:** kinetic (by species), potential (−ε × bound count), and wall heat.
  The first-law residual includes the potential term, as `Simulation.energy()`
  already does for portals.
- **Interventions:** a wall temperature step is a logged command (for E3). ε is
  constant during a run.

## 4. Reference predictions (comparisons only)

Ideal dilute discs; the right chamber is a finite reservoir.

- Free-disc density is the same on both sides: c = c_R.
- Mean bound discs per host: ⟨n⟩ = c·a·e^{ε/T}. The distribution is
  Poisson(⟨n⟩) in the large-reservoir limit, and multinomial when the total
  N_D is fixed (used for the comparison).
- Left total: N_L = c·A_L,acc + N_hosts·⟨n⟩, where A_L,acc excludes
  π(R_out + r)² per host and the posts.
- Osmotic pressure: Π = N_hosts·T / A_L + (host–host second virial term,
  treating hosts as hard discs of radius R_out). The membrane force is
  measured as the net x-impulse on the posts per unit time per unit height.
- μ/T = ln(c_R) on the right, and ln(c_free,L) and ln(n_bound/(a·e^{ε/T})) on
  the left. Shown as differences against a stated reference.

Accessible areas a and A_L,acc are computed analytically (circles) with a Monte
Carlo cross-check. They are stored in each run's metadata. Host motion makes
A_L,acc fluctuate when hosts touch each other or the walls; the cross-check
measures its time average.

Default scene, 200 discs and 8 hosts (ideal-point estimates from the formulas
above, ignoring disc–disc exclusion; to be replaced by measurements):

| ε | T | c | ⟨n⟩ per host | bound | N_left : N_right |
|---:|---:|---:|---:|---:|---:|
| 0 | 1 | 27.4 | 1.0 | 8 | 91 : 109 |
| 1.5 | 0.75 | 21.8 | 5.9 | 47 | **113 : 87** |
| 1.5 | 1 | 24.0 | 4.0 | 32 | 104 : 96 |
| 1.5 | 2 | 26.2 | 2.1 | 16 | **95 : 105** |
| 2.5 | 1 | 18.9 | 8.5 | 68 | 125 : 75 |

The ideal Π/T is 2.0 against a disc pressure P/T of about 24. Π is therefore
about 8 % of the force on each face, which needs long averaging. Fewer discs
make it larger (Π/P = n_hosts/c). E4 uses a lower-density variant.

## 5. Numba engine

**Reused unchanged:** `RunConfig`/presets/CLI, the `World` description, ledger
classes, the thermal-wall kernel and the energy-step formula (moved into
numba-callable functions), `LiveInstruments`, replay archive and window, worker
thread, export, and the test and scientific-check harness. Strictness is
unchanged: an overlap beyond tolerance writes a diagnostic checkpoint and
pauses.

**New:** an event-driven core `core/edmd.py` for "circle worlds" (discs,
rotating rough ring hosts, fixed posts, line walls). The current engine re-plans
every candidate each step from Python, which suits rotating triangles. Here
almost every TOI is an exact quadratic, so a proper event calendar is both simpler and
far faster. A thin `EdmdSimulation` adapter exposes the same surface the worker,
instruments and precalculate already call: `advance_to`, `step_collision`,
`snapshot`, `checkpoint`, `restore`, `save_failure`, `time`, `event_count`,
`world`, `events`/`drain_events`, `apply_command` and the diagnostic counters.

Kernel design, one `@njit(cache=True, fastmath=False)` loop:

- **State:** float64 arrays for x, y, θ, vx, vy, ω, t_last (positions and
  angles are updated only when needed), mass, radius and moment of inertia per
  body, a `kind` array (disc/host),
  per-disc membership (host id or −1), and a per-body collision counter to
  invalidate stale events.
- **Calendar:** a binary heap of (t, type, a, b, feature, counter_a,
  counter_b). Stale entries are dropped when popped.
- **Neighbour search:** a uniform grid with cell size ≥ 2r + margin for
  disc–disc. Each host is registered in every cell its bounding square
  covers and re-registered on its own cell-crossing events. Posts and walls are
  binned once. When a host is hit, only the discs in its registered cells are
  re-predicted. That keeps a host collision O(local), not O(N); this is the
  main design point of the kernel.
- **Disc–host TOI:** solve the circle quadratic for the outer arc (R_out + r),
  the inner arc (R_in − r, from inside) and the step (R_s). Accept an arc root
  only if the contact angle lies outside the mouth at θ(t*). If the disc
  enters the wall band [R_in − r, R_out + r] through the mouth, schedule a
  **mouth-window** search for the end caps over the band interval:
  conservative advancement with the bound |v_rel| + |ω|·R_mid, then
  bisection/Newton to tolerance (§2.1). Host–host, host–post and host–wall use
  the same split: outer-circle quadratic plus a validity check, and a window
  search only when a contact point falls inside a mouth.
- **Contact resolution:** smooth reflection for disc–disc, disc–post and walls;
  rough 2×2 reflection j = −2·D⁻¹·g for contacts with a host (scalar form on
  arcs, where D is diagonal).
- **Random numbers:** a buffer of uniforms drawn from the run's NumPy
  `Generator` and passed in. The kernel returns `NEED_RANDOM` when the buffer
  is empty, so runs are bit-reproducible from the seed.
- **Ledger:** Kahan sums for heat per wall and species, step ΔU by direction,
  membrane x-impulse from each side, entries/exits/reflections per host,
  and the time integral of each occupancy (for averages without sampling bias).
- **Return points:** sample tick, intervention time, `max_events`, buffer
  empty, or an error status. Python never runs inside the loop.

**Speed:** no promise until measured. M1 benchmarks the new core against the
current engine on the plain `gas_box` disc preset (same seed and physics), and
reports events per second for N = 200, 1000 and 5000.

## 6. Experiments and gates

| # | Experiment | Measured | Pass condition |
|---|---|---|---|
| E0a | One host, one trapped disc, gas outside, no membrane | escape-time distribution, step reflections | finite mean escape time; with the host frozen (M, I → ∞) the disc never escapes (shows why mobility matters) |
| E0d | Closed box, hosts and discs, p_rough = 1 and 0.3 | T_trans and T_rot for hosts, T_trans and T_spin for discs | all four agree within error for both p_rough; only the relaxation time differs |
| E0b | Closed box, hosts only, no membrane, ε = 0 and ε > 0 | ⟨n⟩ per host, occupancy histogram, host vs disc temperature | multinomial reference within error; equipartition between species |
| E0c | ε = 0 at T = 0.75 and T = 2 | positional statistics | identical within error |
| E1 | **Membrane**, start with all discs on the right | c_free,L(t), c_R(t), N_L, N_R | free densities converge; N_L ≠ N_R as in §4 |
| E2 | Isotherm: vary N_D at fixed T | ⟨n⟩ vs c_R | the line c·a·e^{ε/T} in the dilute range; departure at high occupancy quantified |
| E3 | **Temperature steps** 0.75 → 2 → 0.75, ε = 1.5 and ε = 0 | N_L − N_R vs time | sign flips with ε = 1.5; no change with ε = 0 |
| E4 | Osmotic force on the fixed membrane, low-density variant, several host counts | Π vs n_hosts·T, Π vs bound count | slope ≈ 1 with the virial term; no dependence on ε at fixed host count |
| E5 | Narrow vs wide mouth, same a and ε | relaxation time after a T step, equilibrium ⟨n⟩ | same equilibrium; different τ |

Each experiment uses several seeds and block averages. Nothing is called
"equilibrium" until the occupancy autocorrelation time is measured and the
run is many times longer.

**Later (needs a moving body under constant force):** the **osmotic piston**.
The membrane is mounted on a guide with a constant opposing force F. It settles
where Π = F/H, so the left chamber visibly swells when hosts are added or
when cooling pulls discs in. Disc–post TOI on a uniformly accelerated
membrane is a quartic over a bounded interval. That is root isolation, not
closed form, so it comes after the fixed-membrane results.

## 7. OsmosisView and instruments

- **Scene:** box, membrane posts, rings drawn with their real mouths and step
  circle (faint), bound discs filled in a second colour, and a brief flash on a
  step crossing (in or out). Optional host trails.
- **Plots:**
  - Free densities (c_free,L and c_R) vs time, with the reference dashed.
  - N_L and N_R vs time, with temperature-step markers.
  - Occupancy histogram per host against the Poisson/multinomial reference.
  - μ_L/T and μ_R/T with block error bars.
  - Membrane force (net) against n_hosts·T.
- **Right column:** T, ε, hosts, bound fraction, Π, ledger residual, and
  "equilibrium reference" labels as in the brief.
- **Replay format v3:** body kind, membership per frame, step events, and the
  osmosis `i_*` columns. Carnot v1/v2 archives are untouched.

## 8. Milestones

| M | Deliverable | Gate |
|---|---|---|
| 1 | `edmd` core for discs + line walls (thermal and specular), adapter, `gas_box` on it | event-for-event parity with the current engine on small N; energy residual ~1e-12; benchmark recorded |
| 2 | Fixed circular posts, grid, membrane geometry checks | analytic TOI tests; discs pass, hosts are blocked in every orientation |
| 3 | Rotating rough ring hosts (arcs, mouth-window caps) without the well | analytic cases incl. a disc entering through a turning mouth; energy, momentum and angular-momentum conservation per contact (~1e-13); E0c, E0d |
| 4 | Binding step with recoil, membership, ledger | climb/reflect/threshold tests; first-law closure with ΔU; E0a, E0b |
| 5 | Presets `osmosis`, `osmosis_hard`; host sizing rule; OsmosisView | E1, E2 |
| 6 | Temperature-step command; osmotic force instrument | E3, E4 |
| 7 | Mouth-width variants; replay v3; video | E5; recorded replay of E3 |
| 8 | Optional: osmotic piston; second inner ring | own gates |

**Progress.** M1 done 2026-10-06 (ADR 0004). Kernel-only throughput is
0.21–0.28 M events/s, nearly flat from N = 200 to 5000; that is 600× the
reference engine at N = 200 and 2700× at N = 1000. Sampling every 0.05 s
with event records costs about 30–40 %. Details are in
`docs/results/edmd_benchmark_2026-10-06.json`. The heap keeps every
prediction, so there is room to optimise (e.g. one event per body) if the
osmosis scenes need it.

M2 done 2026-10-06: fixed circular posts in the kernel and
`experiments/membrane.py` (even gaps, pass/block/mouth checks). Small discs
cross a membrane and spread out evenly; discs larger than the gap stay on
their side over 200 s while hitting the posts. Throughput is now about
0.19 M events/s (ADR 0004).

M3 done 2026-10-06: rotating rough ring hosts (ADR 0004). One change from
§2.1: hosts meet walls, posts and other hosts as their full outer circle
(a "lid" that only discs pass), because a ring's end could otherwise hook
into another ring's mouth. Gate results: per-contact conservation of energy,
momentum and angular momentum; equipartition of all four modes within 3 %
at roughness 1 and 0.3; discs enter and leave cavities; host registration
agrees with a single-cell grid.

M4 done 2026-10-06: binding wells with recoil (ADR 0004). Gate results:
- **Step law:** entering, leaving and refused exits match the law exactly,
  conserving energy, momentum and angular momentum.
- **First law:** closure, potential energy included, at about 1e-11 over
  tens of thousands of crossings.
- **E0a:** a frozen smooth host traps a disc indefinitely, while a mobile
  rough host in a thermal box releases it.
- **E0c:** hard-only occupancy is 0.83–0.87 per well at T = 0.4–2,
  against 0.85 ideal, so independent of T.
- **E0b:** with ε = 1.5, occupancy is 3–7 % below the ideal-point
  reference, with sub-binomial variance. Both come from crowding of
  finite discs in a small well; that is the E2 lesson.

M5 done 2026-10-06: presets `osmosis` and `osmosis_hard`
(`experiments/osmosis.py`, edmd engine by default) with the sizing rule.
More than about 10 hosts make the well smaller than two disc diameters,
which is refused. Also added: `measurements/osmosis.py`, the osmosis panel
in the lab window (`docs/osmosis-lab-1480x900.png`), and removal of the
three stub presets. Gate results:
- **E1:** μ/T on the right equals μ/T of the free left discs within
  0.01–0.05. The hard-only side counts and bound count match the reference
  within 2 %, with fewer discs on the host side. Binding pulls about 9
  more discs onto it.
- **E2:** ⟨n⟩/c_R falls with occupancy (0.138, 0.119, 0.100 at 100, 200,
  400 discs). Extrapolated to empty wells it gives 0.156, against the
  ideal a·e^{ε/T} = 0.166.
- **First look at E4:** the osmotic pressure averaged over 2000 s is 3.5–4.4
  with or without binding, between the ideal n·T = 3.1 and the hard-disc
  virial 5.0.
- **Deferred to M7:** replays of osmosis runs show the rings but not yet
  the membership colours or the osmosis plots (replay v3).

M6 done 2026-10-06. Wall-temperature changes: a `set_wall_temperature`
command (GUI "Wall T"), a schedule (`--temperature-steps 400:2,800:0.75`),
each recorded as a `temperature_change` intervention and restored with
checkpoints.

**E3:** T = 0.75 → 2 → 0.75, 200 discs, two seeds, 300 s windows.

| | excess on host side | bound discs |
|---|---|---|
| ε = 1.5 | +16 to +19 → −4 to −6 → +14 to +18 | 34 → 15 → 34 |
| ε = 0 | −10 to −14 throughout | 8–9 |

**E4:** total net force on the membrane per height, 100 discs, T = 1,
3000 s, block errors ±0.15–0.28.

| hosts | discs | ε | Π | ideal n·T/A |
|---:|---:|---:|---:|---:|
| 8 | 100 | 0 | 3.62, 3.79 | 3.10 |
| 8 | 100 | 1.5 | 3.85, 3.73 | 3.10 |
| 8 | 50 / 200 | 0 | 3.89 / 3.88 | 3.10 |
| 4 | 100 | 0 | 2.21 | 1.88 |
| 2 | 100 | 0 | 1.46 | 1.28 |

Π scales with T and does not depend on the number of discs or on binding.
It sits 15–25 % above the ideal estimate because the hosts exclude each
other. Change from §6: the plan's bulk hard-disc virial reference does not
fit two to eight large hosts in a box and was dropped. The force is also
split by contacting body, but that split is not a pressure of its own:
discs push hosts onto the membrane (fewer fit between a host and the
posts), so the host part grows with disc number (4.08 → 4.46 for 50 → 200
discs) while the total stays.

## 9. Risks

- **Π is small next to the disc pressure.** Mitigated by the low-density E4
  variant and long numba runs; error bars come from blocks.
- **Pockets become crowded at large ε.** Disc–disc exclusion inside a small
  cavity makes the ideal ⟨n⟩ an overestimate. Keep ⟨n⟩ ≲ 5 for the reference
  comparisons, and show the departure (E2) as a lesson.
- **Escape too slow at low T.** The rate falls as e^{−ε/T}, times the mouth
  geometry factor. E0a measures it first; ε and w are tuned before any membrane
  run.
- **Hosts jam against the membrane.** This is harmless physically, but it
  changes A_L,acc. φ = 0.2 leaves room; the accessible-area cross-check reports
  it.
- **Rough contacts are a modelling choice**, like the smooth default
  elsewhere. They are labelled in the UI and metadata, and p_rough = 0 recovers
  smooth rings (with almost no rotational coupling) as a comparison.
- **The mouth-window search is the one non-analytic TOI.** It inherits the
  strict policy: an unresolved bracket pauses the run with a checkpoint.
  Tests include a grazing pass along a cap and a mouth rotating onto a disc
  that sits in the band.

## 10. Decisions (2026-10-06)

1. Hosts **rotate**, coupled by perfectly rough contacts (§2.1).
2. Default binding energy **ε = 1.5**.
3. The stub presets `labyrinth`, `labyrinth_energetic` and `selective_membrane`
   are **removed**, together with `experiments/labyrinth.py`,
   `experiments/permeability.py` and their tests (done in M5, when the
   `osmosis` presets replace them).
