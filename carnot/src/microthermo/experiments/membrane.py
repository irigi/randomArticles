"""Porous membrane: a column of fixed circular posts with checked clearances."""
from __future__ import annotations

from dataclasses import dataclass
import math

from ..core.boundaries import CircularPost


@dataclass(frozen=True)
class Membrane:
    """Posts at x, evenly spaced so every gap (also at the walls) equals ``gap``.

    For circles the narrowest passage between two posts, or between a post
    and a flat wall, lies on the line of closest approach. So a disc of
    radius r passes exactly when 2r < gap.
    """

    x: float
    post_radius: float
    gap: float
    posts: tuple[CircularPost, ...]

    def passes(self, radius: float) -> bool:
        return 2*radius < self.gap


def build_membrane(x: float, y0: float, y1: float, post_radius: float,
                   target_gap: float, passer_radius: float,
                   blocked_radius: float | None = None,
                   mouth_width: float | None = None,
                   clearance: float | None = None) -> Membrane:
    """Choose the post count whose even gap is closest to ``target_gap``.

    Checks, each raising ValueError:
    - discs of ``passer_radius`` pass with ``clearance`` to spare
      (default: a quarter of their radius);
    - bodies of ``blocked_radius`` (e.g. a host's outer radius) cannot pass,
      with the same clearance;
    - a post cannot enter a host mouth of chord ``mouth_width``.
    """
    height = y1 - y0
    if not (height > 0 and post_radius > 0 and target_gap > 0 and passer_radius > 0):
        raise ValueError("membrane needs positive height, post radius, gap and disc radius")
    margin = .25*passer_radius if clearance is None else clearance

    def gap_for(count: int) -> float:
        return (height - 2*count*post_radius)/(count + 1)

    count = max(1, round((height - target_gap)/(2*post_radius + target_gap)))
    candidates = [c for c in (count - 1, count, count + 1) if c >= 1 and gap_for(c) > 0]
    if not candidates:
        raise ValueError("posts do not fit in the membrane height")
    count = min(candidates, key=lambda c: abs(gap_for(c) - target_gap))
    gap = gap_for(count)
    if gap < 2*passer_radius + margin:
        raise ValueError(f"gap {gap:.4g} does not pass discs of radius {passer_radius:.4g} "
                         f"with clearance {margin:.4g}")
    if blocked_radius is not None and gap > 2*blocked_radius - margin:
        raise ValueError(f"gap {gap:.4g} would let bodies of radius {blocked_radius:.4g} pass")
    if mouth_width is not None and 2*post_radius <= mouth_width + margin:
        raise ValueError(f"post diameter {2*post_radius:.4g} could enter a mouth of "
                         f"width {mouth_width:.4g}")
    pitch = 2*post_radius + gap
    posts = tuple(CircularPost((x, y0 + gap + post_radius + k*pitch), post_radius,
                               name=f"membrane_post_{k}")
                  for k in range(count))
    assert math.isclose(posts[-1].center[1] + post_radius + gap, y1, abs_tol=1e-12)
    return Membrane(x, post_radius, gap, posts)
