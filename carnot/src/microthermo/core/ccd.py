from __future__ import annotations

import math
import numpy as np

from .events import Contact, TOIResult, TOIStatus
from .geometry import convex_witnesses, disc_polygon_separation, world_polygon
from .numeric import numba_available, polygon_pair_toi, polygon_witness_at
from .state import BodyArrays, Shape, Tolerances


def disc_disc_toi(pa: np.ndarray, va: np.ndarray, ra: float,
                  pb: np.ndarray, vb: np.ndarray, rb: float, horizon: float,
                  tol: Tolerances = Tolerances()) -> TOIResult:
    dp, dv = pb - pa, vb - va
    r = ra + rb
    c = float(dp @ dp - r*r)
    b = float(dp @ dv)
    if c <= tol.geometry:
        if b >= -tol.velocity:
            return TOIResult(TOIStatus.NO_COLLISION, reason="touching and separating")
        t = 0.0
    else:
        a = float(dv @ dv)
        if a <= tol.velocity**2 or b >= 0.0:
            return TOIResult(TOIStatus.NO_COLLISION)
        disc = b*b - a*c
        if disc < -tol.geometry:
            return TOIResult(TOIStatus.NO_COLLISION)
        disc = max(0.0, disc)
        # cancellation-safe smaller root
        q = -b + math.sqrt(disc)
        t = c / q if q > 0.0 else (-b - math.sqrt(disc)) / a
        if t < -tol.time or t > horizon + tol.time:
            return TOIResult(TOIStatus.NO_COLLISION)
        t = max(0.0, t)
    ca, cb = pa + va*t, pb + vb*t
    d = cb - ca
    norm = float(np.linalg.norm(d))
    if norm <= tol.geometry:
        return TOIResult(TOIStatus.INDETERMINATE, reason="coincident disc centers")
    n = d / norm
    p = ca + ra*n
    return TOIResult(TOIStatus.COLLISION, t, tol.time,
                     Contact(-1, -1, tuple(p), tuple(n)))


def disc_wall_toi(position: np.ndarray, velocity: np.ndarray, radius: float,
                  wall_point: np.ndarray, inward_normal: np.ndarray,
                  horizon: float, wall_velocity: np.ndarray | None = None,
                  tol: Tolerances = Tolerances()) -> TOIResult:
    n = np.asarray(inward_normal, float)
    n /= np.linalg.norm(n)
    wv = np.zeros(2) if wall_velocity is None else np.asarray(wall_velocity, float)
    gap = float((position - wall_point) @ n - radius)
    closing = float((velocity - wv) @ n)
    if gap <= tol.geometry:
        if closing >= -tol.velocity:
            return TOIResult(TOIStatus.NO_COLLISION, reason="touching and separating")
        t = 0.0
    elif closing >= -tol.velocity:
        return TOIResult(TOIStatus.NO_COLLISION)
    else:
        t = -gap / closing
        if t > horizon + tol.time:
            return TOIResult(TOIStatus.NO_COLLISION)
    center = position + velocity*t
    point = center - radius*n
    return TOIResult(TOIStatus.COLLISION, max(t, 0.0), tol.time,
                     Contact(-1, None, tuple(point), tuple(-n)))


def _shape_at(state: BodyArrays, i: int, t: float) -> tuple[str, object]:
    p = state.pos[i] + state.vel[i]*t
    if state.shape[i] == Shape.DISC:
        return "disc", (p, state.radius[i])
    return "poly", world_polygon(state.polygons[i], p, state.angle[i] + state.omega[i]*t)


def _separation(state: BodyArrays, a: int, b: int, t: float,
                numeric_backend: str = "python"):
    if state.shape[a] != Shape.DISC and state.shape[b] != Shape.DISC and (
            numeric_backend == "numba" or
            (numeric_backend == "auto" and numba_available())):
        return polygon_witness_at(
            state.polygons[a], state.polygons[b], state.pos[a], state.pos[b],
            state.vel[a], state.vel[b], state.angle[a], state.angle[b],
            state.omega[a], state.omega[b], t)
    ka, ga = _shape_at(state, a, t)
    kb, gb = _shape_at(state, b, t)
    if ka == kb == "poly":
        sep, n, wa, wb, fa, fb = convex_witnesses(ga, gb)
        return sep, n, .5*(wa+wb), fa, fb
    if ka == "disc" and kb == "poly":
        sep, n, point, _ = disc_polygon_separation(ga[0], ga[1], gb)
        return sep, n, point, -1, -1
    if ka == "poly" and kb == "disc":
        sep, n, point, _ = disc_polygon_separation(gb[0], gb[1], ga)
        return sep, -n, point, -1, -1
    raise AssertionError


def conservative_toi(state: BodyArrays, a: int, b: int, horizon: float,
                     tol: Tolerances = Tolerances(), max_iter: int = 1024,
                     numeric_backend: str = "python") -> TOIResult:
    """Conservative advancement for pairs containing a rotating convex body."""
    if state.shape[a] != Shape.DISC and state.shape[b] != Shape.DISC and (
            numeric_backend == "numba" or
            (numeric_backend == "auto" and numba_available())):
        status, t, error, nx, ny, px, py, fa, fb, reason = polygon_pair_toi(
            state, a, b, horizon, tol, max_iter)
        if status == 1:
            return TOIResult(TOIStatus.COLLISION, t, error,
                             Contact(a, b, (px, py), (nx, ny), fa, fb))
        if status == 0:
            return TOIResult(TOIStatus.NO_COLLISION)
        reasons = ("", "stationary overlap", "initial polygon overlap",
                   "touching contact cannot be isolated",
                   "new contact during touching separation",
                   "touching contact did not separate",
                   "conservative advancement stalled", "iteration limit")
        return TOIResult(TOIStatus.INDETERMINATE, reason=reasons[reason])
    speed_bound = (np.linalg.norm(state.vel[b] - state.vel[a]) +
                   abs(state.omega[a])*state.radius[a] + abs(state.omega[b])*state.radius[b])
    if speed_bound <= tol.velocity:
        sep, *_ = _separation(state, a, b, 0.0, numeric_backend)
        return TOIResult(TOIStatus.NO_COLLISION if sep > tol.geometry else TOIStatus.INDETERMINATE,
                         reason="stationary overlap" if sep <= tol.geometry else "")
    t, prev_t = 0.0, 0.0
    local_time_tol = min(tol.time, tol.geometry/speed_bound)
    for _ in range(max_iter):
        sep, n, point, fa, fb = _separation(state, a, b, t, numeric_backend)
        if sep <= tol.geometry:
            if t <= local_time_tol:
                if sep < -tol.geometry:
                    return TOIResult(TOIStatus.INDETERMINATE, reason="initial polygon overlap")
                ra, rb = point-state.pos[a], point-state.pos[b]
                va = state.vel[a] + state.omega[a]*np.array([-ra[1], ra[0]])
                vb = state.vel[b] + state.omega[b]*np.array([-rb[1], rb[0]])
                if float((vb-va)@n) >= -tol.velocity:
                    # An impact may separate at one feature while rotation
                    # brings a different feature into contact later in the
                    # same interval. Re-query just beyond the current contact.
                    t = min(horizon, max(16*local_time_tol, 16*tol.geometry/speed_bound))
                    if t <= local_time_tol:
                        return TOIResult(TOIStatus.INDETERMINATE, reason="touching contact cannot be isolated")
                    next_sep = _separation(state, a, b, t, numeric_backend)[0]
                    if next_sep < -tol.geometry:
                        return TOIResult(TOIStatus.INDETERMINATE, reason="new contact during touching separation")
                    if next_sep <= tol.geometry:
                        return TOIResult(TOIStatus.INDETERMINATE, reason="touching contact did not separate")
                    prev_t = t
                    continue
            lo, hi = prev_t, t
            for _ in range(60):
                if hi - lo <= local_time_tol:
                    break
                mid = 0.5*(lo + hi)
                if _separation(state, a, b, mid, numeric_backend)[0] <= tol.geometry:
                    hi = mid
                else:
                    lo = mid
            sep, n, point, fa, fb = _separation(state, a, b, hi, numeric_backend)
            return TOIResult(TOIStatus.COLLISION, hi, hi-lo,
                             Contact(a, b, tuple(point), tuple(n), fa, fb))
        prev_t = t
        step = 0.8 * sep / speed_bound
        if step <= local_time_tol:
            return TOIResult(TOIStatus.INDETERMINATE, reason="conservative advancement stalled")
        t += step
        if t > horizon:
            return TOIResult(TOIStatus.NO_COLLISION)
    return TOIResult(TOIStatus.INDETERMINATE, reason="iteration limit")


def body_pair_toi(state: BodyArrays, a: int, b: int, horizon: float,
                  tol: Tolerances = Tolerances(), numeric_backend: str = "python") -> TOIResult:
    if state.shape[a] == Shape.DISC and state.shape[b] == Shape.DISC:
        r = disc_disc_toi(state.pos[a], state.vel[a], state.radius[a],
                          state.pos[b], state.vel[b], state.radius[b], horizon, tol)
        if r.contact:
            c = Contact(a, b, r.contact.point, r.contact.normal)
            return TOIResult(r.status, r.time, r.error, c, r.reason)
        return r
    return conservative_toi(state, a, b, horizon, tol, numeric_backend=numeric_backend)
