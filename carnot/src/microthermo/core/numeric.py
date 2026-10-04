"""Optional Numba kernels for conservative collision queries."""
from __future__ import annotations

import math
import numpy as np

try:
    from numba import njit
except ImportError:  # The reference installation has no compiler dependency.
    njit = None


def _reach_mask_python(pos, vel, radius, omega, pairs, horizon):
    result = np.empty(len(pairs), dtype=np.bool_)
    for k in range(len(pairs)):
        a, b = pairs[k]
        dx = pos[b, 0] - pos[a, 0]
        dy = pos[b, 1] - pos[a, 1]
        vx = vel[b, 0] - vel[a, 0]
        vy = vel[b, 1] - vel[a, 1]
        reach = (radius[a] + radius[b] + horizon *
                 ((vx * vx + vy * vy) ** 0.5 +
                  abs(omega[a]) * radius[a] + abs(omega[b]) * radius[b]))
        result[k] = dx * dx + dy * dy <= reach * reach
    return result


if njit is not None:
    _reach_mask_numba = njit(cache=True, fastmath=False)(_reach_mask_python)
else:
    _reach_mask_numba = None


def reach_mask(pos, vel, radius, omega, pairs, horizon, backend="auto"):
    """Filter grid candidates with the same conservative reach as the reference."""
    if backend not in ("auto", "python", "numba"):
        raise ValueError("numeric backend must be auto, python, or numba")
    if backend == "numba" and _reach_mask_numba is None:
        raise RuntimeError("Numba backend requires installation of microthermo[accel]")
    fn = _reach_mask_numba if backend != "python" and _reach_mask_numba else _reach_mask_python
    return fn(pos, vel, radius, omega, pairs, horizon)


def _swept_circle_mask_python(pos, vel, radius, pairs, horizon, geom_tol):
    """Exact closest-center approach of linearly moving circumcircles."""
    result = np.empty(len(pairs), dtype=np.bool_)
    for k in range(len(pairs)):
        a, b = pairs[k]
        dx = pos[b, 0] - pos[a, 0]
        dy = pos[b, 1] - pos[a, 1]
        vx = vel[b, 0] - vel[a, 0]
        vy = vel[b, 1] - vel[a, 1]
        speed_sq = vx * vx + vy * vy
        t = min(horizon, max(0.0, -(dx * vx + dy * vy) / speed_sq)) if speed_sq > 0 else 0.0
        x = dx + vx * t
        y = dy + vy * t
        # The geometry margin also protects near-tangent floating-point cases.
        reach = radius[a] + radius[b] + 2.0 * geom_tol
        result[k] = x * x + y * y <= reach * reach
    return result


_swept_circle_mask_numba = (njit(cache=True, fastmath=False)(_swept_circle_mask_python)
                            if njit is not None else None)


def swept_circle_mask(pos, vel, radius, pairs, horizon, geom_tol,
                      backend="auto"):
    """Keep pairs whose moving circumcircles may meet in the interval."""
    if backend not in ("auto", "python", "numba"):
        raise ValueError("numeric backend must be auto, python, or numba")
    if backend == "numba" and _swept_circle_mask_numba is None:
        raise RuntimeError("Numba backend requires installation of microthermo[accel]")
    fn = (_swept_circle_mask_numba if backend != "python" and
          _swept_circle_mask_numba is not None else _swept_circle_mask_python)
    return fn(pos, vel, radius, pairs, horizon, geom_tol)


def _wall_reach_mask_python(pos, vel, radius, wall_points, wall_normals,
                            wall_speeds, segment_mask, horizon, geom_tol):
    result = np.empty((len(radius), len(wall_speeds)), dtype=np.bool_)
    for a in range(len(radius)):
        for wi in range(len(wall_speeds)):
            if segment_mask[wi]:
                result[a, wi] = True
                continue
            nx, ny = wall_normals[wi, 0], wall_normals[wi, 1]
            gap = ((pos[a, 0] - wall_points[wi, 0]) * nx +
                   (pos[a, 1] - wall_points[wi, 1]) * ny - radius[a])
            speed = abs(vel[a, 0] * nx + vel[a, 1] * ny)
            # A wider margin makes this mask a conservative superset of the
            # scalar wall bound in Simulation._wall_toi.
            result[a, wi] = gap <= horizon * (speed + wall_speeds[wi]) + 2 * geom_tol
    return result


_wall_reach_mask_numba = (njit(cache=True, fastmath=False)(_wall_reach_mask_python)
                          if njit is not None else None)


def wall_reach_mask(pos, vel, radius, wall_points, wall_normals, wall_speeds,
                    segment_mask, horizon, geom_tol, backend="auto"):
    """Conservative batch prefilter for particle-wall interval queries."""
    if backend not in ("auto", "python", "numba"):
        raise ValueError("numeric backend must be auto, python, or numba")
    if backend == "numba" and _wall_reach_mask_numba is None:
        raise RuntimeError("Numba backend requires installation of microthermo[accel]")
    fn = (_wall_reach_mask_numba if backend != "python" and
          _wall_reach_mask_numba is not None else _wall_reach_mask_python)
    return fn(pos, vel, radius, wall_points, wall_normals, wall_speeds,
              segment_mask, horizon, geom_tol)


def numba_available() -> bool:
    return _reach_mask_numba is not None


def _polygon_witness_raw(a, b):
    """Scalar SAT and closest-feature witness for convex polygon pairs."""
    na, nb = len(a), len(b)
    cax = cay = cbx = cby = 0.0
    for i in range(na):
        cax += a[i, 0]
        cay += a[i, 1]
    for i in range(nb):
        cbx += b[i, 0]
        cby += b[i, 1]
    cax /= na
    cay /= na
    cbx /= nb
    cby /= nb
    separation = -math.inf
    axis_x, axis_y = 1.0, 0.0
    for group in range(2):
        poly = a if group == 0 else b
        n = len(poly)
        for i in range(n):
            nxt = (i + 1) % n
            ex = poly[nxt, 0] - poly[i, 0]
            ey = poly[nxt, 1] - poly[i, 1]
            norm = math.sqrt(ex * ex + ey * ey)
            if norm == 0.0:
                continue
            nx, ny = -ey / norm, ex / norm
            if (cbx - cax) * nx + (cby - cay) * ny < 0.0:
                nx, ny = -nx, -ny
            amin = math.inf
            amax = -math.inf
            bmin = math.inf
            for j in range(na):
                projection = a[j, 0] * nx + a[j, 1] * ny
                amin = min(amin, projection)
                amax = max(amax, projection)
            for j in range(nb):
                projection = b[j, 0] * nx + b[j, 1] * ny
                bmin = min(bmin, projection)
            sep = bmin - amax
            if sep > separation:
                separation, axis_x, axis_y = sep, nx, ny
    best = math.inf
    wax = way = wbx = wby = 0.0
    fa = fb = -1
    for i in range(na):
        vx, vy = a[i, 0], a[i, 1]
        for j in range(nb):
            sx, sy = b[j, 0], b[j, 1]
            ex = b[(j + 1) % nb, 0] - sx
            ey = b[(j + 1) % nb, 1] - sy
            dd = ex * ex + ey * ey
            u = max(0.0, min(1.0, ((vx - sx) * ex + (vy - sy) * ey) / dd)) if dd else 0.0
            qx, qy = sx + u * ex, sy + u * ey
            dx, dy = qx - vx, qy - vy
            d2 = dx * dx + dy * dy
            if d2 < best:
                best, wax, way, wbx, wby, fa = d2, vx, vy, qx, qy, 2 * i
                fb = 2 * j if u == 0.0 else (2 * ((j + 1) % nb) if u == 1.0 else 2 * j + 1)
    for j in range(nb):
        vx, vy = b[j, 0], b[j, 1]
        for i in range(na):
            sx, sy = a[i, 0], a[i, 1]
            ex = a[(i + 1) % na, 0] - sx
            ey = a[(i + 1) % na, 1] - sy
            dd = ex * ex + ey * ey
            u = max(0.0, min(1.0, ((vx - sx) * ex + (vy - sy) * ey) / dd)) if dd else 0.0
            qx, qy = sx + u * ex, sy + u * ey
            dx, dy = qx - vx, qy - vy
            d2 = dx * dx + dy * dy
            if d2 < best:
                best, wax, way, wbx, wby, fb = d2, qx, qy, vx, vy, 2 * j
                fa = 2 * i if u == 0.0 else (2 * ((i + 1) % na) if u == 1.0 else 2 * i + 1)
    length = math.sqrt(best)
    if separation > 0.0 and fb % 2 == 1:
        j = fb // 2
        ex = b[(j + 1) % nb, 0] - b[j, 0]
        ey = b[(j + 1) % nb, 1] - b[j, 1]
        area = 0.0
        for i in range(nb):
            area += b[i, 0] * b[(i + 1) % nb, 1] - b[i, 1] * b[(i + 1) % nb, 0]
        orientation = 1.0 if area > 0.0 else -1.0
        norm = math.sqrt(ex * ex + ey * ey)
        nx, ny = orientation * -ey / norm, orientation * ex / norm
    elif separation > 0.0 and fa % 2 == 1:
        i = fa // 2
        ex = a[(i + 1) % na, 0] - a[i, 0]
        ey = a[(i + 1) % na, 1] - a[i, 1]
        area = 0.0
        for j in range(na):
            area += a[j, 0] * a[(j + 1) % na, 1] - a[j, 1] * a[(j + 1) % na, 0]
        orientation = 1.0 if area > 0.0 else -1.0
        norm = math.sqrt(ex * ex + ey * ey)
        nx, ny = orientation * ey / norm, orientation * -ex / norm
    elif separation > 0.0 and length > 0.0:
        nx, ny = (wbx - wax) / length, (wby - way) / length
    else:
        nx, ny = axis_x, axis_y
    return separation, nx, ny, wax, way, wbx, wby, fa, fb


if njit is not None:
    _polygon_witness_numba = njit(cache=True, fastmath=False)(_polygon_witness_raw)
else:
    _polygon_witness_numba = None


def polygon_witnesses(a, b, backend="auto"):
    if backend == "numba" and _polygon_witness_numba is None:
        raise RuntimeError("Numba backend requires installation of microthermo[accel]")
    fn = _polygon_witness_numba if backend != "python" and _polygon_witness_numba else _polygon_witness_raw
    sep, nx, ny, ax, ay, bx, by, fa, fb = fn(a, b)
    return sep, np.array([nx, ny]), np.array([ax, ay]), np.array([bx, by]), fa, fb


if njit is not None:
    @njit(cache=True, fastmath=False)
    def _polygon_sat_gap_numba(a, b):
        """The SAT gap from the witness kernel, without closest features."""
        na, nb = len(a), len(b)
        cax = cay = cbx = cby = 0.0
        for i in range(na):
            cax += a[i, 0]
            cay += a[i, 1]
        for i in range(nb):
            cbx += b[i, 0]
            cby += b[i, 1]
        cax /= na
        cay /= na
        cbx /= nb
        cby /= nb
        separation = -math.inf
        for group in range(2):
            poly = a if group == 0 else b
            n = len(poly)
            for i in range(n):
                nxt = (i + 1) % n
                ex = poly[nxt, 0] - poly[i, 0]
                ey = poly[nxt, 1] - poly[i, 1]
                norm = math.sqrt(ex * ex + ey * ey)
                if norm == 0.0:
                    continue
                nx, ny = -ey / norm, ex / norm
                if (cbx - cax) * nx + (cby - cay) * ny < 0.0:
                    nx, ny = -nx, -ny
                amin = math.inf
                amax = -math.inf
                bmin = math.inf
                for j in range(na):
                    projection = a[j, 0] * nx + a[j, 1] * ny
                    amin = min(amin, projection)
                    amax = max(amax, projection)
                for j in range(nb):
                    projection = b[j, 0] * nx + b[j, 1] * ny
                    bmin = min(bmin, projection)
                sep = bmin - amax
                if sep > separation:
                    separation = sep
        return separation

    @njit(cache=True, fastmath=False)
    def _polygon_gap_at_numba(local_a, local_b, pos_a, pos_b, vel_a, vel_b,
                              angle_a, angle_b, omega_a, omega_b, t):
        aa, ab = angle_a + omega_a * t, angle_b + omega_b * t
        ca, sa = math.cos(aa), math.sin(aa)
        cb, sb = math.cos(ab), math.sin(ab)
        a = np.empty((len(local_a), 2), dtype=np.float64)
        b = np.empty((len(local_b), 2), dtype=np.float64)
        for i in range(len(local_a)):
            a[i, 0] = local_a[i, 0] * ca - local_a[i, 1] * sa + pos_a[0] + vel_a[0] * t
            a[i, 1] = local_a[i, 0] * sa + local_a[i, 1] * ca + pos_a[1] + vel_a[1] * t
        for i in range(len(local_b)):
            b[i, 0] = local_b[i, 0] * cb - local_b[i, 1] * sb + pos_b[0] + vel_b[0] * t
            b[i, 1] = local_b[i, 0] * sb + local_b[i, 1] * cb + pos_b[1] + vel_b[1] * t
        return _polygon_sat_gap_numba(a, b)

    @njit(cache=True, fastmath=False)
    def _polygon_witness_at_numba(local_a, local_b, pos_a, pos_b, vel_a, vel_b,
                                  angle_a, angle_b, omega_a, omega_b, t):
        aa = angle_a + omega_a * t
        ab = angle_b + omega_b * t
        ca, sa = math.cos(aa), math.sin(aa)
        cb, sb = math.cos(ab), math.sin(ab)
        a = np.empty((len(local_a), 2), dtype=np.float64)
        b = np.empty((len(local_b), 2), dtype=np.float64)
        for i in range(len(local_a)):
            a[i, 0] = local_a[i, 0] * ca - local_a[i, 1] * sa + pos_a[0] + vel_a[0] * t
            a[i, 1] = local_a[i, 0] * sa + local_a[i, 1] * ca + pos_a[1] + vel_a[1] * t
        for i in range(len(local_b)):
            b[i, 0] = local_b[i, 0] * cb - local_b[i, 1] * sb + pos_b[0] + vel_b[0] * t
            b[i, 1] = local_b[i, 0] * sb + local_b[i, 1] * cb + pos_b[1] + vel_b[1] * t
        return _polygon_witness_numba(a, b)
else:
    _polygon_witness_at_numba = None


def polygon_witness_at(local_a, local_b, pos_a, pos_b, vel_a, vel_b,
                       angle_a, angle_b, omega_a, omega_b, t):
    sep, nx, ny, ax, ay, bx, by, fa, fb = _polygon_witness_at_numba(
        local_a, local_b, pos_a, pos_b, vel_a, vel_b,
        angle_a, angle_b, omega_a, omega_b, t)
    return sep, np.array([nx, ny]), np.array([(ax+bx)/2, (ay+by)/2]), fa, fb


if njit is not None:
    @njit(cache=True, fastmath=False)
    def _polygon_pair_vertices_into_numba(local_a, local_b, pos_a, pos_b,
                                          vel_a, vel_b, angle_a, angle_b,
                                          omega_a, omega_b, t, a, b):
        aa, ab = angle_a + omega_a * t, angle_b + omega_b * t
        ca, sa = math.cos(aa), math.sin(aa)
        cb, sb = math.cos(ab), math.sin(ab)
        for i in range(len(local_a)):
            a[i, 0] = local_a[i, 0] * ca - local_a[i, 1] * sa + pos_a[0] + vel_a[0] * t
            a[i, 1] = local_a[i, 0] * sa + local_a[i, 1] * ca + pos_a[1] + vel_a[1] * t
        for i in range(len(local_b)):
            b[i, 0] = local_b[i, 0] * cb - local_b[i, 1] * sb + pos_b[0] + vel_b[0] * t
            b[i, 1] = local_b[i, 0] * sb + local_b[i, 1] * cb + pos_b[1] + vel_b[1] * t

    @njit(cache=True, fastmath=False)
    def _fixed_polygon_swept_separated_numba(a, b, dvx, dvy, horizon, geom_tol):
        for poly in (a, b):
            for i in range(len(poly)):
                j = (i + 1) % len(poly)
                nx = -(poly[j, 1] - poly[i, 1])
                ny = poly[j, 0] - poly[i, 0]
                norm = math.sqrt(nx * nx + ny * ny)
                if norm == 0.0:
                    continue
                amin, amax = math.inf, -math.inf
                bmin, bmax = math.inf, -math.inf
                for k in range(len(a)):
                    p = a[k, 0] * nx + a[k, 1] * ny
                    amin, amax = min(amin, p), max(amax, p)
                for k in range(len(b)):
                    p = b[k, 0] * nx + b[k, 1] * ny
                    bmin, bmax = min(bmin, p), max(bmax, p)
                travel = (dvx * nx + dvy * ny) * horizon
                if (bmin + min(0.0, travel) - amax > geom_tol * norm or
                        amin - bmax - max(0.0, travel) > geom_tol * norm):
                    return True
        return False

    @njit(cache=True, fastmath=False)
    def _polygon_pair_toi_numba(local_a, local_b, pos_a, pos_b, vel_a, vel_b,
                                angle_a, angle_b, omega_a, omega_b,
                                radius_a, radius_b, horizon, geom_tol,
                                time_tol, velocity_tol, max_iter, gap_only,
                                a, b):
        dvx = vel_b[0] - vel_a[0]
        dvy = vel_b[1] - vel_a[1]
        speed_bound = (math.sqrt(dvx * dvx + dvy * dvy) +
                       abs(omega_a) * radius_a + abs(omega_b) * radius_b)
        if omega_a == 0.0 and omega_b == 0.0:
            _polygon_pair_vertices_into_numba(local_a, local_b, pos_a, pos_b,
                                              vel_a, vel_b, angle_a, angle_b,
                                              omega_a, omega_b, 0.0, a, b)
            if _fixed_polygon_swept_separated_numba(
                    a, b, dvx, dvy, horizon, geom_tol):
                return 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1, -1, 0
        def gap_at(t):
            _polygon_pair_vertices_into_numba(local_a, local_b, pos_a, pos_b,
                                              vel_a, vel_b, angle_a, angle_b,
                                              omega_a, omega_b, t, a, b)
            return _polygon_sat_gap_numba(a, b)
        def witness_at(t):
            _polygon_pair_vertices_into_numba(local_a, local_b, pos_a, pos_b,
                                              vel_a, vel_b, angle_a, angle_b,
                                              omega_a, omega_b, t, a, b)
            return _polygon_witness_numba(a, b)
        if speed_bound <= velocity_tol:
            if gap_only:
                sep = gap_at(0.0)
            else:
                sep = witness_at(0.0)[0]
            return (0 if sep > geom_tol else 2, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, -1, -1,
                    0 if sep > geom_tol else 1)
        t = 0.0
        prev_t = 0.0
        local_time_tol = min(time_tol, geom_tol / speed_bound)
        for _ in range(max_iter):
            if gap_only:
                sep = gap_at(t)
            if not gap_only or sep <= geom_tol + 1e-14:
                sep, nx, ny, ax, ay, bx, by, fa, fb = witness_at(t)
            if sep <= geom_tol:
                if t <= local_time_tol:
                    if sep < -geom_tol:
                        return 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1, -1, 2
                    px, py = .5 * (ax + bx), .5 * (ay + by)
                    rax, ray = px - pos_a[0], py - pos_a[1]
                    rbx, rby = px - pos_b[0], py - pos_b[1]
                    vax = vel_a[0] - omega_a * ray
                    vay = vel_a[1] + omega_a * rax
                    vbx = vel_b[0] - omega_b * rby
                    vby = vel_b[1] + omega_b * rbx
                    if (vbx - vax) * nx + (vby - vay) * ny >= -velocity_tol:
                        t = min(horizon, max(16 * local_time_tol,
                                             16 * geom_tol / speed_bound))
                        if t <= local_time_tol:
                            return 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1, -1, 3
                        if gap_only:
                            next_sep = gap_at(t)
                        else:
                            next_sep = witness_at(t)[0]
                        if next_sep < -geom_tol:
                            return 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1, -1, 4
                        if next_sep <= geom_tol:
                            return 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1, -1, 5
                        prev_t = t
                        continue
                lo, hi = prev_t, t
                for _ in range(60):
                    if hi - lo <= local_time_tol:
                        break
                    mid = .5 * (lo + hi)
                    if gap_only:
                        mid_sep = gap_at(mid)
                    else:
                        mid_sep = witness_at(mid)[0]
                    if mid_sep <= geom_tol:
                        hi = mid
                    else:
                        lo = mid
                sep, nx, ny, ax, ay, bx, by, fa, fb = witness_at(hi)
                return (1, hi, hi - lo, nx, ny,
                        .5 * (ax + bx), .5 * (ay + by), fa, fb, 0)
            prev_t = t
            step = .8 * sep / speed_bound
            if step <= local_time_tol:
                return 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1, -1, 6
            t += step
            if t > horizon:
                return 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1, -1, 0
        return 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1, -1, 7
else:
    _polygon_pair_toi_numba = None


def polygon_pair_toi(state, a, b, horizon, tol, max_iter):
    """Return a raw pair TOI; status 0=no hit, 1=hit, 2=indeterminate."""
    return _polygon_pair_toi_numba(
        state.polygons[a], state.polygons[b], state.pos[a], state.pos[b],
        state.vel[a], state.vel[b], state.angle[a], state.angle[b],
        state.omega[a], state.omega[b], state.radius[a], state.radius[b],
        horizon, tol.geometry, tol.time, tol.velocity, max_iter, False,
        np.empty((len(state.polygons[a]), 2), dtype=np.float64),
        np.empty((len(state.polygons[b]), 2), dtype=np.float64))


if njit is not None:
    @njit(cache=True, fastmath=False)
    def _polygon_pairs_toi_numba(local, pos, vel, angle, omega, radius,
                                 pairs, horizon, geom_tol, time_tol,
                                 velocity_tol, max_iter, gap_only):
        result = np.empty((len(pairs), 10), dtype=np.float64)
        scratch_a = np.empty((local.shape[1], 2), dtype=np.float64)
        scratch_b = np.empty((local.shape[1], 2), dtype=np.float64)
        for k in range(len(pairs)):
            a, b = pairs[k, 0], pairs[k, 1]
            raw = _polygon_pair_toi_numba(
                local[a], local[b], pos[a], pos[b], vel[a], vel[b],
                angle[a], angle[b], omega[a], omega[b], radius[a], radius[b],
                horizon, geom_tol, time_tol, velocity_tol, max_iter, gap_only,
                scratch_a, scratch_b)
            result[k, 0] = raw[0]
            result[k, 1] = raw[1]
            result[k, 2] = raw[2]
            result[k, 3] = raw[3]
            result[k, 4] = raw[4]
            result[k, 5] = raw[5]
            result[k, 6] = raw[6]
            result[k, 7] = raw[7]
            result[k, 8] = raw[8]
            result[k, 9] = raw[9]
        return result
else:
    _polygon_pairs_toi_numba = None


def polygon_pairs_toi(state, pairs, horizon, tol, max_iter, gap_only=True):
    """Batch polygon pair CCD while preserving the input pair order."""
    if _polygon_pairs_toi_numba is None:
        raise RuntimeError("Numba backend requires installation of microthermo[accel]")
    local = (state.packed_triangles if state.packed_triangles is not None else
             np.stack([state.polygons[i] for i in range(state.n)]))
    return _polygon_pairs_toi_numba(
        local, state.pos, state.vel, state.angle, state.omega, state.radius,
        pairs, horizon, tol.geometry, tol.time, tol.velocity, max_iter,
        gap_only)


if njit is not None:
    @njit(cache=True, fastmath=False)
    def _triangle_pair_gaps_numba(local, pos, vel, angle, omega, radius, pairs):
        """Current SAT gaps for sorted triangle candidates."""
        result = np.empty(len(pairs), dtype=np.float64)
        for k in range(len(pairs)):
            a, b = pairs[k, 0], pairs[k, 1]
            dx = pos[b, 0] - pos[a, 0]
            dy = pos[b, 1] - pos[a, 1]
            if math.sqrt(dx * dx + dy * dy) > radius[a] + radius[b]:
                result[k] = math.inf
            else:
                result[k] = _polygon_gap_at_numba(
                    local[a], local[b], pos[a], pos[b], vel[a], vel[b],
                    angle[a], angle[b], omega[a], omega[b], 0.0)
        return result
else:
    _triangle_pair_gaps_numba = None


def triangle_pair_gaps(state, pairs):
    """Batch overlap gaps for triangle-only current-state diagnostics."""
    if _triangle_pair_gaps_numba is None:
        raise RuntimeError("Numba backend requires installation of microthermo[accel]")
    local = (state.packed_triangles if state.packed_triangles is not None else
             np.stack([state.polygons[i] for i in range(state.n)]))
    return _triangle_pair_gaps_numba(local, state.pos, state.vel, state.angle,
                                     state.omega, state.radius, pairs)


if njit is not None:
    @njit(cache=True, fastmath=False)
    def _disc_pairs_toi_numba(pos, vel, radius, pairs, horizon,
                              geom_tol, time_tol, velocity_tol):
        result = np.zeros((len(pairs), 8), dtype=np.float64)
        for k in range(len(pairs)):
            a, b = pairs[k, 0], pairs[k, 1]
            dx, dy = pos[b, 0] - pos[a, 0], pos[b, 1] - pos[a, 1]
            vx, vy = vel[b, 0] - vel[a, 0], vel[b, 1] - vel[a, 1]
            r = radius[a] + radius[b]
            c = dx * dx + dy * dy - r * r
            approach = dx * vx + dy * vy
            if c <= geom_tol:
                if approach >= -velocity_tol:
                    continue
                t = 0.0
            else:
                speed_sq = vx * vx + vy * vy
                if speed_sq <= velocity_tol * velocity_tol or approach >= 0.0:
                    continue
                discriminant = approach * approach - speed_sq * c
                if discriminant < -geom_tol:
                    continue
                discriminant = max(0.0, discriminant)
                q = -approach + math.sqrt(discriminant)
                t = c / q if q > 0.0 else (-approach - math.sqrt(discriminant)) / speed_sq
                if t < -time_tol or t > horizon + time_tol:
                    continue
                t = max(0.0, t)
            cax, cay = pos[a, 0] + vel[a, 0] * t, pos[a, 1] + vel[a, 1] * t
            cbx, cby = pos[b, 0] + vel[b, 0] * t, pos[b, 1] + vel[b, 1] * t
            nx, ny = cbx - cax, cby - cay
            norm = math.sqrt(nx * nx + ny * ny)
            if norm <= geom_tol:
                result[k, 0] = 2.0
                continue
            nx, ny = nx / norm, ny / norm
            result[k, 0] = 1.0
            result[k, 1] = t
            result[k, 2] = time_tol
            result[k, 3] = cax + radius[a] * nx
            result[k, 4] = cay + radius[a] * ny
            result[k, 5] = nx
            result[k, 6] = ny
        return result
else:
    _disc_pairs_toi_numba = None


def disc_pairs_toi(state, pairs, horizon, tol):
    """Batch disc CCD; status 0=no hit, 1=hit, 2=coincident centers."""
    if _disc_pairs_toi_numba is None:
        raise RuntimeError("Numba backend requires installation of microthermo[accel]")
    return _disc_pairs_toi_numba(state.pos, state.vel, state.radius, pairs,
                                 horizon, tol.geometry, tol.time, tol.velocity)


if njit is not None:
    @njit(cache=True, fastmath=False)
    def _controlled_cam_profile(phi, areas, boundaries, height):
        phase = phi % (2.0 * math.pi)
        sector = 3
        for i in range(3):
            if phase < boundaries[i + 1]:
                sector = i
                break
        width = boundaries[sector + 1] - boundaries[sector]
        u = (phase - boundaries[sector]) / width
        u = min(1.0, max(0.0, u))
        delta = areas[(sector + 1) % 4] - areas[sector]
        smooth = u**3 * (10.0 + u * (-15.0 + 6.0 * u))
        derivative = 30.0 * u * u * (1.0 - u) * (1.0 - u)
        return ((areas[sector] + delta * smooth) / height,
                delta * derivative / (width * height))

    @njit(cache=True, fastmath=False)
    def _controlled_cam_disc_toi_numba(pos, vel, radius, possible,
                                       phase0, shaft_omega, areas, boundaries,
                                       height, wall_speed, horizon,
                                       geom_tol, time_tol, velocity_tol):
        # Rows: status, time, error, point x, point y, failure reason.
        result = np.zeros((len(radius), 6), dtype=np.float64)
        for a in range(len(radius)):
            if not possible[a]:
                continue
            bound = (math.sqrt(vel[a, 0]**2 + vel[a, 1]**2) + wall_speed)
            if bound <= velocity_tol:
                continue
            local_time_tol = min(time_tol, geom_tol / bound)
            t = previous = 0.0
            for _ in range(1024):
                wall_x, slope = _controlled_cam_profile(
                    phase0 + shaft_omega * t, areas, boundaries, height)
                gap = wall_x - (pos[a, 0] + vel[a, 0] * t) - radius[a]
                relative = slope * shaft_omega - vel[a, 0]
                if gap <= geom_tol:
                    if t == 0.0 and relative >= -velocity_tol:
                        break
                    lo, hi = previous, t
                    for _ in range(60):
                        if hi - lo <= local_time_tol:
                            break
                        mid = 0.5 * (lo + hi)
                        mid_x = _controlled_cam_profile(
                            phase0 + shaft_omega * mid, areas, boundaries, height)[0]
                        mid_gap = mid_x - (pos[a, 0] + vel[a, 0] * mid) - radius[a]
                        if mid_gap <= geom_tol:
                            hi = mid
                        else:
                            lo = mid
                    result[a, 0] = 1.0
                    result[a, 1] = hi
                    result[a, 2] = hi - lo
                    result[a, 3] = pos[a, 0] + vel[a, 0] * hi + radius[a]
                    result[a, 4] = pos[a, 1] + vel[a, 1] * hi
                    break
                previous = t
                step = 0.8 * gap / bound
                if step <= local_time_tol:
                    result[a, 0] = 2.0
                    result[a, 5] = 1.0
                    break
                t += step
                if t > horizon:
                    break
            else:
                result[a, 0] = 2.0
                result[a, 5] = 2.0
        return result
else:
    _controlled_cam_disc_toi_numba = None


def controlled_cam_disc_toi(state, possible, mechanism, wall_speed, horizon, tol):
    """Batch controlled-shaft cam piston contacts for discs."""
    if _controlled_cam_disc_toi_numba is None:
        raise RuntimeError("Numba backend requires installation of microthermo[accel]")
    return _controlled_cam_disc_toi_numba(
        state.pos, state.vel, state.radius, possible,
        mechanism.shaft.phi, mechanism.shaft.prescribed_omega,
        np.asarray(mechanism.cam.areas), mechanism.cam.boundaries,
        mechanism.cam.height, wall_speed, horizon,
        tol.geometry, tol.time, tol.velocity)


if njit is not None:
    @njit(cache=True, fastmath=False)
    def _cam_triangle_gap(local, pos, vel, angle, omega,
                          phase0, shaft_omega, areas, boundaries, height, t):
        theta = angle + omega * t
        c, s = math.cos(theta), math.sin(theta)
        cx, cy = pos[0] + vel[0] * t, pos[1] + vel[1] * t
        best_x = -math.inf
        px = py = 0.0
        feature = -1
        for i in range(len(local)):
            x = local[i, 0] * c - local[i, 1] * s + cx
            y = local[i, 0] * s + local[i, 1] * c + cy
            if x > best_x:
                best_x, px, py, feature = x, x, y, i
        wall_x, slope = _controlled_cam_profile(
            phase0 + shaft_omega * t, areas, boundaries, height)
        return wall_x - best_x, px, py, feature, slope * shaft_omega

    @njit(cache=True, fastmath=False)
    def _controlled_cam_triangle_toi_numba(local, pos, vel, angle, omega,
                                           radius, possible, phase0,
                                           shaft_omega, areas, boundaries,
                                           height, wall_speed, horizon,
                                           geom_tol, time_tol, velocity_tol):
        # Rows: status, time, error, point x, point y, vertex, reason.
        result = np.zeros((len(radius), 7), dtype=np.float64)
        for a in range(len(radius)):
            if not possible[a]:
                continue
            bound = (math.sqrt(vel[a, 0]**2 + vel[a, 1]**2) +
                     abs(omega[a]) * radius[a] + wall_speed)
            if bound <= velocity_tol:
                continue
            local_time_tol = min(time_tol, geom_tol / bound)
            t = previous = 0.0
            finished = False
            for _ in range(1024):
                gap, px, py, feature, wall_vx = _cam_triangle_gap(
                    local[a], pos[a], vel[a], angle[a], omega[a],
                    phase0, shaft_omega, areas, boundaries, height, t)
                if gap <= geom_tol:
                    if gap < -geom_tol and t == 0.0:
                        result[a, 0], result[a, 6] = 2.0, 1.0
                        finished = True
                        break
                    cy = pos[a, 1] + vel[a, 1] * t
                    relative = wall_vx - vel[a, 0] + omega[a] * (py - cy)
                    if t <= local_time_tol and relative >= -velocity_tol:
                        t = min(horizon, max(16 * local_time_tol,
                                             16 * geom_tol / bound))
                        if t <= local_time_tol:
                            result[a, 0], result[a, 6] = 2.0, 2.0
                            finished = True
                            break
                        next_gap = _cam_triangle_gap(
                            local[a], pos[a], vel[a], angle[a], omega[a],
                            phase0, shaft_omega, areas, boundaries, height, t)[0]
                        if next_gap < -geom_tol:
                            result[a, 0], result[a, 6] = 2.0, 3.0
                            finished = True
                            break
                        if next_gap <= geom_tol:
                            result[a, 0], result[a, 6] = 2.0, 4.0
                            finished = True
                            break
                        previous = t
                        continue
                    lo, hi = previous, t
                    for _ in range(60):
                        if hi - lo <= local_time_tol:
                            break
                        mid = .5 * (lo + hi)
                        mid_gap = _cam_triangle_gap(
                            local[a], pos[a], vel[a], angle[a], omega[a],
                            phase0, shaft_omega, areas, boundaries, height, mid)[0]
                        if mid_gap <= geom_tol:
                            hi = mid
                        else:
                            lo = mid
                    _, px, py, feature, _ = _cam_triangle_gap(
                        local[a], pos[a], vel[a], angle[a], omega[a],
                        phase0, shaft_omega, areas, boundaries, height, hi)
                    result[a, 0] = 1.0
                    result[a, 1] = hi
                    result[a, 2] = hi - lo
                    result[a, 3] = px
                    result[a, 4] = py
                    result[a, 5] = feature
                    finished = True
                    break
                previous = t
                step = .8 * gap / bound
                if step <= local_time_tol:
                    result[a, 0], result[a, 6] = 2.0, 5.0
                    finished = True
                    break
                t += step
                if t > horizon:
                    finished = True
                    break
            if not finished:
                result[a, 0], result[a, 6] = 2.0, 6.0
        return result
else:
    _controlled_cam_triangle_toi_numba = None


def controlled_cam_triangle_toi(state, possible, mechanism, wall_speed,
                                horizon, tol):
    """Batch controlled-shaft cam piston contacts for triangles."""
    if _controlled_cam_triangle_toi_numba is None:
        raise RuntimeError("Numba backend requires installation of microthermo[accel]")
    local = (state.packed_triangles if state.packed_triangles is not None else
             np.stack([state.polygons[i] for i in range(state.n)]))
    return _controlled_cam_triangle_toi_numba(
        local, state.pos, state.vel, state.angle, state.omega, state.radius,
        possible, mechanism.shaft.phi, mechanism.shaft.prescribed_omega,
        np.asarray(mechanism.cam.areas), mechanism.cam.boundaries,
        mechanism.cam.height, wall_speed, horizon,
        tol.geometry, tol.time, tol.velocity)


if njit is not None:
    @njit(cache=True, fastmath=False)
    def _fixed_wall_gap(local, pos, vel, angle, omega, nx, ny,
                        wall_projection, wall_normal_speed, t):
        theta = angle + omega * t
        c, s = math.cos(theta), math.sin(theta)
        cx = pos[0] + vel[0] * t
        cy = pos[1] + vel[1] * t
        best = math.inf
        px = py = 0.0
        feature = -1
        for i in range(len(local)):
            x = local[i, 0] * c - local[i, 1] * s + cx
            y = local[i, 0] * s + local[i, 1] * c + cy
            projection = x * nx + y * ny
            if projection < best:
                best, px, py, feature = projection, x, y, i
        return best - (wall_projection + wall_normal_speed * t), px, py, feature

    @njit(cache=True, fastmath=False)
    def _fixed_wall_toi_numba(local, pos, vel, angle, omega, radius,
                              inward, wall_start, wall_velocity, speed,
                              horizon, geom_tol, time_tol, velocity_tol):
        if speed <= velocity_tol:
            return 0, 0.0, 0.0, 0.0, 0.0, -1, 0
        nx, ny = inward[0], inward[1]
        wall_projection = wall_start[0] * nx + wall_start[1] * ny
        wall_normal_speed = wall_velocity[0] * nx + wall_velocity[1] * ny
        local_time_tol = min(time_tol, geom_tol / speed)
        t = previous = 0.0
        for _ in range(1024):
            gap, px, py, feature = _fixed_wall_gap(
                local, pos, vel, angle, omega, nx, ny,
                wall_projection, wall_normal_speed, t)
            if gap <= geom_tol:
                if gap < -geom_tol and t == 0.0:
                    return 2, 0.0, 0.0, 0.0, 0.0, -1, 1
                arm_x = px - pos[0] - vel[0] * t
                arm_y = py - pos[1] - vel[1] * t
                relative_x = vel[0] - wall_velocity[0] - omega * arm_y
                relative_y = vel[1] - wall_velocity[1] + omega * arm_x
                if t <= local_time_tol and relative_x * nx + relative_y * ny >= -velocity_tol:
                    t = min(horizon, max(16 * local_time_tol, 16 * geom_tol / speed))
                    if t <= local_time_tol:
                        return 2, 0.0, 0.0, 0.0, 0.0, -1, 2
                    next_gap = _fixed_wall_gap(
                        local, pos, vel, angle, omega, nx, ny,
                        wall_projection, wall_normal_speed, t)[0]
                    if next_gap < -geom_tol:
                        return 2, 0.0, 0.0, 0.0, 0.0, -1, 3
                    if next_gap <= geom_tol:
                        return 2, 0.0, 0.0, 0.0, 0.0, -1, 4
                    previous = t
                    continue
                lo, hi = previous, t
                for _ in range(60):
                    if hi - lo <= local_time_tol:
                        break
                    mid = .5 * (lo + hi)
                    if _fixed_wall_gap(local, pos, vel, angle, omega,
                                       nx, ny, wall_projection,
                                       wall_normal_speed, mid)[0] <= geom_tol:
                        hi = mid
                    else:
                        lo = mid
                _, px, py, feature = _fixed_wall_gap(
                    local, pos, vel, angle, omega, nx, ny,
                    wall_projection, wall_normal_speed, hi)
                return 1, hi, hi - lo, px, py, feature, 0
            previous = t
            step = .8 * gap / speed
            if step <= local_time_tol:
                return 2, 0.0, 0.0, 0.0, 0.0, -1, 5
            t += step
            if t > horizon:
                return 0, 0.0, 0.0, 0.0, 0.0, -1, 0
        return 2, 0.0, 0.0, 0.0, 0.0, -1, 6
else:
    _fixed_wall_toi_numba = None


def fixed_wall_toi(state, body, wall, wall_start, speed, horizon, tol):
    """Raw polygon TOI for a wall with constant linear velocity."""
    return _fixed_wall_toi_numba(
        state.polygons[body], state.pos[body], state.vel[body],
        state.angle[body], state.omega[body], state.radius[body],
        wall.inward_normal, wall_start, wall.velocity, speed, horizon,
        tol.geometry, tol.time, tol.velocity)


if njit is not None:
    @njit(cache=True, fastmath=False)
    def _fixed_triangle_walls_toi_numba(local, pos, vel, angle, omega, radius,
                                        possible, starts, normals, wall_velocities,
                                        wall_speeds, horizon, geom_tol,
                                        time_tol, velocity_tol):
        result = np.zeros((len(radius), len(starts), 7), dtype=np.float64)
        for body in range(len(radius)):
            vx, vy = vel[body, 0], vel[body, 1]
            for wi in range(len(starts)):
                if not possible[body, wi]:
                    continue
                nx, ny = normals[wi, 0], normals[wi, 1]
                gap_bound = ((pos[body, 0] - starts[wi, 0]) * nx +
                             (pos[body, 1] - starts[wi, 1]) * ny - radius[body])
                body_speed = abs(vx * nx + vy * ny)
                if gap_bound > horizon * (body_speed + wall_speeds[wi]) + geom_tol:
                    continue
                speed = (math.sqrt(vx * vx + vy * vy) +
                         abs(omega[body]) * radius[body] + wall_speeds[wi])
                raw = _fixed_wall_toi_numba(
                    local[body], pos[body], vel[body], angle[body], omega[body],
                    radius[body], normals[wi], starts[wi], wall_velocities[wi],
                    speed, horizon, geom_tol, time_tol, velocity_tol)
                result[body, wi, 0] = raw[0]
                result[body, wi, 1] = raw[1]
                result[body, wi, 2] = raw[2]
                result[body, wi, 3] = raw[3]
                result[body, wi, 4] = raw[4]
                result[body, wi, 5] = raw[5]
                result[body, wi, 6] = raw[6]
        return result
else:
    _fixed_triangle_walls_toi_numba = None


def fixed_triangle_walls_toi(state, possible, starts, normals, velocities,
                             speeds, horizon, tol):
    """Batch triangle TOI against constant-velocity infinite-line walls."""
    if _fixed_triangle_walls_toi_numba is None:
        raise RuntimeError("Numba backend requires installation of microthermo[accel]")
    local = (state.packed_triangles if state.packed_triangles is not None else
             np.stack([state.polygons[i] for i in range(state.n)]))
    return _fixed_triangle_walls_toi_numba(
        local, state.pos, state.vel, state.angle, state.omega, state.radius,
        possible, starts, normals, velocities, speeds, horizon,
        tol.geometry, tol.time, tol.velocity)


if njit is not None:
    @njit(cache=True, fastmath=False)
    def _fixed_disc_walls_toi_numba(pos, vel, radius, possible, starts,
                                    normals, wall_velocities, wall_speeds,
                                    horizon, geom_tol, time_tol, velocity_tol):
        """Reference disc advancement for several constant-velocity walls."""
        result = np.zeros((len(radius), len(starts), 6), dtype=np.float64)
        for body in range(len(radius)):
            for wi in range(len(starts)):
                if not possible[body, wi]:
                    continue
                nx, ny = normals[wi, 0], normals[wi, 1]
                vx, vy = vel[body, 0], vel[body, 1]
                bound = math.sqrt(vx * vx + vy * vy) + wall_speeds[wi]
                if bound <= velocity_tol:
                    continue
                local_time_tol = min(time_tol, geom_tol / bound)
                wall_projection = starts[wi, 0] * nx + starts[wi, 1] * ny
                wall_normal_speed = wall_velocities[wi, 0] * nx + wall_velocities[wi, 1] * ny
                t = previous = 0.0
                finished = False
                for _ in range(1024):
                    cx = pos[body, 0] + vx * t
                    cy = pos[body, 1] + vy * t
                    gap = (cx * nx + cy * ny - wall_projection -
                           wall_normal_speed * t - radius[body])
                    rel = ((vx - wall_velocities[wi, 0]) * nx +
                           (vy - wall_velocities[wi, 1]) * ny)
                    if gap <= geom_tol:
                        if t == 0.0 and rel >= -velocity_tol:
                            finished = True
                            break
                        lo, hi = previous, t
                        for _ in range(60):
                            if hi - lo <= local_time_tol:
                                break
                            mid = .5 * (lo + hi)
                            mx = pos[body, 0] + vx * mid
                            my = pos[body, 1] + vy * mid
                            mid_gap = (mx * nx + my * ny - wall_projection -
                                       wall_normal_speed * mid - radius[body])
                            if mid_gap <= geom_tol:
                                hi = mid
                            else:
                                lo = mid
                        cx = pos[body, 0] + vx * hi
                        cy = pos[body, 1] + vy * hi
                        result[body, wi, 0] = 1.0
                        result[body, wi, 1] = hi
                        result[body, wi, 2] = hi - lo
                        result[body, wi, 3] = cx - radius[body] * nx
                        result[body, wi, 4] = cy - radius[body] * ny
                        finished = True
                        break
                    previous = t
                    step = .8 * gap / bound
                    if step <= local_time_tol:
                        result[body, wi, 0] = 2.0
                        result[body, wi, 5] = 1.0
                        finished = True
                        break
                    t += step
                    if t > horizon:
                        finished = True
                        break
                if not finished:
                    result[body, wi, 0] = 2.0
                    result[body, wi, 5] = 2.0
        return result
else:
    _fixed_disc_walls_toi_numba = None


def fixed_disc_walls_toi(state, possible, starts, normals, velocities,
                         speeds, horizon, tol):
    """Batch disc TOI against constant-velocity infinite-line walls."""
    if _fixed_disc_walls_toi_numba is None:
        raise RuntimeError("Numba backend requires installation of microthermo[accel]")
    return _fixed_disc_walls_toi_numba(
        state.pos, state.vel, state.radius, possible, starts, normals,
        velocities, speeds, horizon, tol.geometry, tol.time, tol.velocity)


if njit is not None:
    @njit(cache=True, fastmath=False)
    def _wall_lower_bounds_numba(pos, radius, wall_points, wall_normals):
        bounds = np.empty((len(radius), len(wall_points)), dtype=np.float64)
        for i in range(len(radius)):
            for j in range(len(wall_points)):
                dx = pos[i, 0] - wall_points[j, 0]
                dy = pos[i, 1] - wall_points[j, 1]
                bounds[i, j] = (dx * wall_normals[j, 0] +
                                dy * wall_normals[j, 1] - radius[i])
        return bounds
else:
    _wall_lower_bounds_numba = None


def wall_lower_bounds(pos, radius, wall_points, wall_normals):
    """Circumcircle lower bounds for all body and infinite-line wall pairs."""
    return _wall_lower_bounds_numba(pos, radius, wall_points, wall_normals)


if njit is not None:
    @njit(cache=True, fastmath=False)
    def _disc_penetration_numba(pos, radius, pairs, wall_points, wall_normals):
        """Visit sorted disc pairs and infinite walls in body-major order."""
        deepest = 0.0
        kind = body_index = other_index = -1
        pair_index = 0
        for a in range(len(radius)):
            while pair_index < len(pairs) and pairs[pair_index, 0] == a:
                b = pairs[pair_index, 1]
                dx = pos[b, 0] - pos[a, 0]
                dy = pos[b, 1] - pos[a, 1]
                gap = math.sqrt(dx * dx + dy * dy) - radius[a] - radius[b]
                if gap < deepest:
                    deepest, kind, body_index, other_index = gap, 1, a, b
                pair_index += 1
            for wi in range(len(wall_points)):
                dx = pos[a, 0] - wall_points[wi, 0]
                dy = pos[a, 1] - wall_points[wi, 1]
                gap = dx * wall_normals[wi, 0] + dy * wall_normals[wi, 1] - radius[a]
                if gap < deepest:
                    deepest, kind, body_index, other_index = gap, 2, a, wi
        return deepest, kind, body_index, other_index
else:
    _disc_penetration_numba = None


def disc_penetration(pos, radius, pairs, wall_points, wall_normals):
    """Deepest overlap and participants for disc-only infinite-wall worlds."""
    if _disc_penetration_numba is None:
        raise RuntimeError("Numba backend requires installation of microthermo[accel]")
    return _disc_penetration_numba(pos, radius, pairs, wall_points, wall_normals)


if njit is not None:
    @njit(cache=True, fastmath=False)
    def _triangle_penetration_numba(local, pos, vel, angle, omega, radius,
                                    pairs, wall_points, wall_normals):
        """Visit sorted triangle pairs and infinite walls in body-major order."""
        deepest = 0.0
        kind = body_index = other_index = -1
        pair_index = 0
        for a in range(len(radius)):
            while pair_index < len(pairs) and pairs[pair_index, 0] == a:
                b = pairs[pair_index, 1]
                dx = pos[b, 0] - pos[a, 0]
                dy = pos[b, 1] - pos[a, 1]
                if math.sqrt(dx * dx + dy * dy) <= radius[a] + radius[b]:
                    gap = _polygon_gap_at_numba(
                        local[a], local[b], pos[a], pos[b], vel[a], vel[b],
                        angle[a], angle[b], omega[a], omega[b], 0.0)
                    if gap < deepest:
                        deepest, kind, body_index, other_index = gap, 1, a, b
                pair_index += 1
            for wi in range(len(wall_points)):
                nx, ny = wall_normals[wi, 0], wall_normals[wi, 1]
                projection = wall_points[wi, 0] * nx + wall_points[wi, 1] * ny
                gap = _fixed_wall_gap(local[a], pos[a], vel[a], angle[a],
                                      omega[a], nx, ny, projection, 0.0, 0.0)[0]
                if gap < deepest:
                    deepest, kind, body_index, other_index = gap, 2, a, wi
        return deepest, kind, body_index, other_index
else:
    _triangle_penetration_numba = None


def triangle_penetration(state, pairs, wall_points, wall_normals):
    """Deepest overlap and participants for triangle-only infinite-wall worlds."""
    if _triangle_penetration_numba is None:
        raise RuntimeError("Numba backend requires installation of microthermo[accel]")
    local = (state.packed_triangles if state.packed_triangles is not None else
             np.stack([state.polygons[i] for i in range(state.n)]))
    return _triangle_penetration_numba(
        local, state.pos, state.vel, state.angle, state.omega, state.radius,
        pairs, wall_points, wall_normals)


def polygon_wall_gap(state, body, wall_point, inward):
    """Exact polygon support gap at the current state, using the numeric kernel."""
    projection = wall_point[0] * inward[0] + wall_point[1] * inward[1]
    return _fixed_wall_gap(
        state.polygons[body], state.pos[body], state.vel[body],
        state.angle[body], state.omega[body], inward[0], inward[1],
        projection, 0.0, 0.0)[0]
