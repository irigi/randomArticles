"""Conservative candidate enumeration for bounded straight-line motion."""
from __future__ import annotations

from collections import defaultdict
import math
import numpy as np

try:
    from numba import njit
except ImportError:
    njit = None


def swept_pairs(pos: np.ndarray, vel: np.ndarray, radius: np.ndarray,
                horizon: float, cell_size: float | None = None) -> list[tuple[int, int]]:
    """Return every pair whose swept circumcircle boxes can intersect.

    A rotating polygon always remains inside its circumcircle. The box covers
    the complete center trajectory over the horizon, so omitted pairs cannot
    collide, including at a grid boundary.
    """
    n = len(radius)
    if n < 2:
        return []
    size = cell_size or max(4.0 * float(np.max(radius)), 0.05)
    cells: dict[tuple[int, int], list[int]] = defaultdict(list)
    for i in range(n):
        end_x = pos[i, 0] + horizon * vel[i, 0]
        end_y = pos[i, 1] + horizon * vel[i, 1]
        lo_x = math.floor((min(pos[i, 0], end_x) - radius[i]) / size)
        hi_x = math.floor((max(pos[i, 0], end_x) + radius[i]) / size)
        lo_y = math.floor((min(pos[i, 1], end_y) - radius[i]) / size)
        hi_y = math.floor((max(pos[i, 1], end_y) + radius[i]) / size)
        for x in range(lo_x, hi_x + 1):
            for y in range(lo_y, hi_y + 1):
                cells[x, y].append(i)
    pairs: set[tuple[int, int]] = set()
    for members in cells.values():
        for index, a in enumerate(members):
            for b in members[index + 1:]:
                pairs.add((a, b))
    return sorted(pairs)


def all_pairs(n: int) -> list[tuple[int, int]]:
    return [(a, b) for a in range(n) for b in range(a + 1, n)]


def _swept_pairs_sweep_raw(pos, vel, radius, horizon):
    """Enumerate overlapping swept circumcircle boxes by their x intervals."""
    n = len(radius)
    lo_x = np.empty(n)
    hi_x = np.empty(n)
    lo_y = np.empty(n)
    hi_y = np.empty(n)
    for i in range(n):
        ex = pos[i, 0] + horizon * vel[i, 0]
        ey = pos[i, 1] + horizon * vel[i, 1]
        lo_x[i] = min(pos[i, 0], ex) - radius[i]
        hi_x[i] = max(pos[i, 0], ex) + radius[i]
        lo_y[i] = min(pos[i, 1], ey) - radius[i]
        hi_y[i] = max(pos[i, 1], ey) + radius[i]
    order = np.argsort(lo_x)
    count = 0
    for k in range(n):
        a = order[k]
        for m in range(k + 1, n):
            b = order[m]
            if lo_x[b] > hi_x[a]:
                break
            if lo_y[b] <= hi_y[a] and lo_y[a] <= hi_y[b]:
                count += 1
    pairs = np.empty((count, 2), dtype=np.int64)
    index = 0
    for k in range(n):
        a = order[k]
        for m in range(k + 1, n):
            b = order[m]
            if lo_x[b] > hi_x[a]:
                break
            if lo_y[b] <= hi_y[a] and lo_y[a] <= hi_y[b]:
                pairs[index, 0] = min(a, b)
                pairs[index, 1] = max(a, b)
                index += 1
    return pairs


_swept_pairs_sweep_numba = (njit(cache=True, fastmath=False)(_swept_pairs_sweep_raw)
                            if njit is not None else None)


def swept_pairs_sweep_array(pos: np.ndarray, vel: np.ndarray, radius: np.ndarray,
                            horizon: float, backend: str = "auto") -> np.ndarray:
    """Return sorted x/y swept-box candidates as an (n, 2) integer array."""
    if backend not in ("auto", "python", "numba"):
        raise ValueError("numeric backend must be auto, python, or numba")
    if backend == "numba" and _swept_pairs_sweep_numba is None:
        raise RuntimeError("Numba backend requires installation of microthermo[accel]")
    fn = (_swept_pairs_sweep_numba if backend != "python" and
          _swept_pairs_sweep_numba is not None else _swept_pairs_sweep_raw)
    pairs = fn(pos, vel, radius, horizon)
    if len(pairs):
        pairs = pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))]
    return pairs


def swept_pairs_sweep(pos: np.ndarray, vel: np.ndarray, radius: np.ndarray,
                      horizon: float, backend: str = "auto") -> list[tuple[int, int]]:
    """Return x/y swept-box candidates in stable body-index order."""
    pairs = swept_pairs_sweep_array(pos, vel, radius, horizon, backend)
    return [tuple(map(int, pair)) for pair in pairs]
