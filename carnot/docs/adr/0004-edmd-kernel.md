# ADR 0004: compiled event-driven kernel for circle worlds

Status: accepted, 2026-10-06 (osmosis plan, milestone 1).

The reference engine (ADR 0001) re-plans every candidate time of impact from
Python at each bounded-horizon step. That suits rotating triangles, whose TOI
needs conservative advancement. For worlds made of circles every TOI is an
exact quadratic, so a classical event-driven kernel is simpler and much
faster.

`core/edmd.py` is one Numba loop:

- **Calendar:** binary heap of all predicted events, invalidated by per-body
  collision counters. It is compacted by a full replan when it nears
  capacity, and Python doubles the capacity if that is not enough.
- **Lazy positions:** body i is at `pos[i] + vel[i]*(t - tl[i])`. Sampling
  evaluates positions without touching the stored state, so the sample
  cadence cannot change the trajectory (tested bitwise).
- **Neighbour search:** a uniform grid with cells of at least the largest
  diameter. Cell-crossing events update membership; after a crossing only
  the newly adjacent cells are predicted.
- **Random numbers:** Python passes a buffer of uniforms drawn from the run's
  `Generator`. One draw per thermal-wall contact, in event order, so a seed
  consumes the same stream as in the reference engine.
- **Strictness:** an overlap found while predicting, or at a sample, raises
  `NumericalFailure` with a checkpoint, as in the reference engine.

`runner/edmd_simulation.EdmdSimulation` exposes the reference `Simulation`
interface. It is selected with `RunConfig.engine = "edmd"` (`--engine edmd`)
and rejects worlds it does not support.

**Parity.** Disc–disc times agree with the reference engine to about 1e-16.
Wall contact times differ by about 1e-10, because the reference engine
finds wall contacts by conservative advancement to a 1e-10 gap. Chaos then
amplifies the difference, so event sequences agree for the first few dozen
events (tested), after which only statistics can be compared.

**Fixed circular posts** (milestone 2) are specular static circles
(`World.posts`, reference engine refuses them). Each post is listed in every
cell its bounding square, grown by the largest disc radius, overlaps. A disc
tests the posts of its own cell when it enters the cell and after every
velocity change; any disc touching a post has its centre in a listed cell,
so no contact is missed (tested against a single-cell grid). Impulse on each
post is Kahan-summed for the osmotic-force instrument.

**Arrays travel in five tuples** (bodies, grid, static geometry, calendar,
I/O) with the helpers inlined. This costs about 15–20 % against the
milestone-1 kernel's long argument lists (0.19 vs 0.24 M events/s at
N = 1000, same machine), accepted for readability. The larger lever, if
needed, is keeping one pending event per body instead of every prediction.

**Rotating ring hosts** (milestone 3). A host is a *thick arc*: every point
within h/2 of an arc of radius R_mid that leaves out the mouth. The distance
from a disc to the host is its distance to that arc minus h/2. Inside the
arc's span that distance is radial, so outer and inner surface contacts are
quadratics that rotation does not change; only the rounded mouth ends
rotate. A disc can touch them only while crossing the wall band through the
mouth, an interval given by two more quadratics. Only inside it does a
conservative-advancement search run (bound |v_rel| + |omega|·R_mid; an
unresolved search is a strict failure).

- **Rough contacts.** Disc–host, host–host and host–post contacts are
  perfectly rough with probability `World.contact_roughness`: the whole
  relative contact velocity reverses. The 2×2 effective-mass inverse keeps
  energy, momentum and angular momentum (tested per contact to ~1e-13).
  Discs carry spin (I = m r²/2), which only these contacts change.
  Equipartition between disc translation and spin and host translation and
  rotation holds within about 3 % at roughness 1 and 0.3.
- **Mouth lid (modelling choice).** Hosts meet walls, posts and other hosts
  as their full outer circle, as if a lid that only discs can pass covered
  the mouth. With a 3-diameter mouth and a 0.03 wall, one ring's end could
  otherwise slide into another's mouth and hook the two together. The lid
  rules that out, keeps those contacts quadratic, and makes the host–host
  reference the hard-disc virial. Discs see the real mouth.
- **Locality.** Hosts are not in the disc grid. Each host is registered in
  the cells around it (outer circle plus the largest disc radius plus a
  two-cell margin) and re-registers when it has moved by the margin, so a
  host hit re-predicts only nearby discs (tested against a single cell).
- **Throughput.** About 0.15 M events/s for plain discs and 0.13 M with
  six hosts and 200 discs, against 0.19 M at milestone 2 (same session).
  The cost is the larger compiled loop; per-function inlining did not
  recover it. The structural fix, if needed, is one pending event per body.

**Binding wells** (milestone 4). A ring may carry a concentric step circle
(`RingGeometry.well_radius`, `well_depth`). Crossing it is a quadratic event
resolved by the brief's reversible step law (section 8.3) with host recoil:
D = 1/m + 1/M, an uphill step below the 1e-12 threshold of
`resolve_energy_step` reflects. The impulse is radial through both centres,
so it exerts no torque. The kernel tracks well membership, entries, exits and
refused exits, and the potential energy, which the first-law residual
includes. With depth 0 the circle only counts occupancy.

Scope: discs (spinless or spinning) and rotating rough ring hosts with
optional binding wells, in an axis-aligned box of stationary specular or
thermal walls, with fixed posts.
