# ADR 0001: strict NumPy reference backend first

Status: accepted, 2026-09-28.

The first executable backend uses Float64 NumPy state and a deterministic
bounded-horizon scheduler. Analytic disc queries and conservative advancement
for rotating convex polygons return collision, proven miss, or indeterminate.
An indeterminate result pauses the run with a reproducible error.

This keeps collision order and energy accounting inspectable on the available
Python environment, which has NumPy but not Numba or Qt. Numba remains an
optional dependency and should be introduced only after profiling while the
reference functions stay available for cross-checks. Fast-math is excluded
from scientific kernels.

