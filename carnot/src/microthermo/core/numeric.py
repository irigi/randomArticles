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
    def _polygon_pair_toi_numba(local_a, local_b, pos_a, pos_b, vel_a, vel_b,
                                angle_a, angle_b, omega_a, omega_b,
                                radius_a, radius_b, horizon, geom_tol,
                                time_tol, velocity_tol, max_iter):
        dvx = vel_b[0] - vel_a[0]
        dvy = vel_b[1] - vel_a[1]
        speed_bound = (math.sqrt(dvx * dvx + dvy * dvy) +
                       abs(omega_a) * radius_a + abs(omega_b) * radius_b)
        if speed_bound <= velocity_tol:
            sep = _polygon_witness_at_numba(
                local_a, local_b, pos_a, pos_b, vel_a, vel_b,
                angle_a, angle_b, omega_a, omega_b, 0.0)[0]
            return (0 if sep > geom_tol else 2, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, -1, -1,
                    0 if sep > geom_tol else 1)
        t = 0.0
        prev_t = 0.0
        local_time_tol = min(time_tol, geom_tol / speed_bound)
        for _ in range(max_iter):
            sep, nx, ny, ax, ay, bx, by, fa, fb = _polygon_witness_at_numba(
                local_a, local_b, pos_a, pos_b, vel_a, vel_b,
                angle_a, angle_b, omega_a, omega_b, t)
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
                        next_sep = _polygon_witness_at_numba(
                            local_a, local_b, pos_a, pos_b, vel_a, vel_b,
                            angle_a, angle_b, omega_a, omega_b, t)[0]
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
                    mid_sep = _polygon_witness_at_numba(
                        local_a, local_b, pos_a, pos_b, vel_a, vel_b,
                        angle_a, angle_b, omega_a, omega_b, mid)[0]
                    if mid_sep <= geom_tol:
                        hi = mid
                    else:
                        lo = mid
                sep, nx, ny, ax, ay, bx, by, fa, fb = _polygon_witness_at_numba(
                    local_a, local_b, pos_a, pos_b, vel_a, vel_b,
                    angle_a, angle_b, omega_a, omega_b, hi)
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
        horizon, tol.geometry, tol.time, tol.velocity, max_iter)


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


def polygon_wall_gap(state, body, wall_point, inward):
    """Exact polygon support gap at the current state, using the numeric kernel."""
    projection = wall_point[0] * inward[0] + wall_point[1] * inward[1]
    return _fixed_wall_gap(
        state.polygons[body], state.pos[body], state.vel[body],
        state.angle[body], state.omega[body], inward[0], inward[1],
        projection, 0.0, 0.0)[0]
