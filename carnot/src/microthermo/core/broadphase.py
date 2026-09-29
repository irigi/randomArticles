"""Conservative candidate enumeration for bounded straight-line motion."""
from __future__ import annotations

from collections import defaultdict
import math
import numpy as np


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
