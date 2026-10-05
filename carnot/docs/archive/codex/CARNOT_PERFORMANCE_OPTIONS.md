# Carnot performance options for review

Prepared: 2026-09-29. **Decision: A+B selected.** The first grid and Numba
filter slice is implemented; see `CARNOT_GAP_CLOSURE_PLAN.md` for measured
results and remaining work. The baseline below predates particle-radius scaling.

## Measured baseline

These are single unprofiled runs on an Intel Core i7-11850H, CPython 3.12.3,
NumPy 2.4.6, seed 123, `max_horizon=0.05`, and the current Python reference
backend. Initialization is excluded. The short intervals measure solver cost,
not long-run scientific throughput or a guaranteed UI frame rate.

| Preset | Particles | Physical time | Wall time | Events | Physical time / wall time |
| --- | ---: | ---: | ---: | ---: | ---: |
| Carnot discs | 16 | 1.0 | 0.199 s | 48 | 5.02× |
| Carnot discs | 32 | 1.0 | 1.278 s | 136 | 0.78× |
| Carnot discs | 64 | 0.5 | 6.406 s | 240 | 0.078× |
| Carnot triangles | 16 | 0.5 | 1.448 s | 31 | 0.345× |
| Carnot triangles | 32 | 0.5 | 9.046 s | 78 | 0.055× |

A `cProfile` pass on 32 discs assigned 1.23 of 2.13 profiled seconds to
`_earliest` and 0.78 seconds to `_penetration`. Both scan all pairs repeatedly;
`_earliest` also calls every wall query for every particle. On 16 triangles,
polygon pair TOI took 1.69 of 2.01 profiled seconds; its repeated feature
witness queries dominate. These categories contain nested functions and must
not be summed with their children. Profiling adds overhead, so use the first
table for absolute baseline timing.

The cam's `boundaries` property rebuilds a NumPy array repeatedly, and small
2D distance tests create many NumPy temporaries. These are credible local
costs, but the all-pairs scan and triangle TOI are the larger targets.

The former Carnot preset fixed radius at `0.025` and minimum cylinder area at `1.25`.
At 500 discs the particle area alone is about `500*pi*0.025^2/1.25 = 0.785`
of the minimum cylinder area. Larger-particle benchmarks need an explicit
radius or geometry scaling rule; otherwise initialization and the ideal-gas
comparison become physically misleading before CPU speed is addressed.
The chosen rule is `radius(N) = 0.025*sqrt(48/N)`. It preserves the area
fraction at the 48-particle reference: about 3.12% for triangles and 7.54%
for discs in the minimum cylinder. It changes collision rates and still
increases gas pressure with N. Particle placement still checks each new
candidate against previously placed particles.

## Options

| Route | What it changes | Advantages | Costs and risks |
| --- | --- | --- | --- |
| A. Better search, Python reference | Swept spatial grid, deterministic candidate order, cached cam/wall invariants, fewer tiny array allocations. | No new build tool; reduces unnecessary pair work for discs and triangles; benefits every backend. | Safe swept bounds for rotating polygons and moving piston need careful tests; Python polygon TOI still limits throughput. |
| B. A plus Numba kernels | Compile numeric broad phase, penetration, disc TOI, then triangle witness/CCD loops with `@njit(cache=True, fastmath=False)`. Keep event orchestration, RNG, ledger, and GUI in Python. | Fits existing NumPy arrays; incremental optional backend; good match to profiled numerical loops. | Requires numeric kernel interfaces rather than decorating object-heavy `Simulation`; JIT warm-up, version pinning, and Python fallback need maintenance. |
| C. A plus selective Rust kernels | Use PyO3/maturin for batch broad phase and CCD calls over arrays, keeping Python orchestration. | Predictable native build/runtime; can release the interpreter lock around large batches; a route to deeper optimization. | More code and wheel packaging; per-pair Python/Rust calls would erase gains, so a coarse batch API is required. Geometry and event-order parity still need full validation. |
| D. Full Rust simulation core | Port scheduler, mechanism, collisions, ledger, and RNG. | Maximum control over data layout and scheduling. | Highest rewrite and validation cost; changing RNG/floating-point/event order complicates parity. Consider only after measured A–C results show a need. |

The selected A+B path is incremental. Numba kernels now handle swept reach
filtering and the full polygon pair CCD loop, including transforms, witnesses,
touching checks, and bisection. Conservative circle bounds prune most wall
queries; Numba now handles fixed-velocity polygon wall CCD and batched
penetration bounds. Event orchestration, cam piston motion, and the ledger
remain Python.
Consider C only if mature A+B measurements miss the 200/500-particle target.

## Compatibility and environment

- `numba` is declared as an optional `accel` dependency and version 0.67.0 is
  installed in the current `.venv`; `requirements-accel-lock.txt` pins this
  tested pair with NumPy 2.4.6. Numba's current
  published compatibility table lists version 0.67.0 as supporting Python
  3.10–3.14 and NumPy 2.0–2.5, including this environment. Before using B,
  verify a clean install of the tested pair on the target machine.
  Numba's nopython mode is required for the intended speedup; `fastmath` is disabled.
- `rustc 1.96.0` and `cargo 1.96.0` are already on this machine. `maturin` is
  absent. If C is chosen, I can set up a separate optional Rust crate and a
  local development install. The expected starting commands are:

  ```bash
  source .venv/bin/activate
  python -m pip install maturin
  maturin develop --release --manifest-path rust-kernels/Cargo.toml
  ```

  The crate, binding API, build configuration, and installation check would be
  created and tested as part of that option. No Rust environment changes are
  needed to compare the options now.

## Acceptance contract shared by all routes

1. Keep the reference Python backend available. Run fixed-seed parity tests
   for event type, chronological order, participants, contact geometry,
   numerical-health counters, heat/work, RNG state, and final snapshots.
2. Cover discs, rotating triangles, thermal and moving walls, exact branch
   boundaries, near-simultaneous contacts, checkpoint/resume, and irregular
   sampling. Declare tolerances before comparing. A speedup cannot hide CCD
   failures or penetration.
3. Benchmark cold startup separately from warm throughput. Publish hardware,
   Python/NumPy/backend versions, seed, particle radius/density, particle count,
   physical duration, wall time, event count, memory, and GUI achieved rate.
4. Keep candidate enumeration deterministic and conservative. A spatial grid
   may reject only pairs proven unable to meet during the search horizon.
   Event invalidation may be added later only after scheduler parity holds.

## Official references

- [Numba version support](https://numba.readthedocs.io/en/stable/user/installing.html#version-support-information)
- [Numba performance guidance](https://numba.readthedocs.io/en/stable/user/performance-tips.html)
- [Numba JIT options](https://numba.readthedocs.io/en/stable/reference/jit-compilation.html)
- [PyO3 getting started](https://pyo3.rs/main/getting-started)
- [PyO3 performance and interpreter detachment](https://pyo3.rs/main/performance.html)
- [Maturin local development](https://www.maturin.rs/local_development)
