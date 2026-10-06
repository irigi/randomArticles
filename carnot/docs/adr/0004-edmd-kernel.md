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

Scope at milestone 1: smooth spinless discs in an axis-aligned box of
stationary specular or thermal walls. Fixed circular posts, rotating rough
ring hosts and the binding step follow in milestones 2–4 of
`docs/PLAN_OSMOSIS.md`.
