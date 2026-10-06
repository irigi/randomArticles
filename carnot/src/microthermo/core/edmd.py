"""Event-driven kernel for smooth discs among fixed walls and circular posts.

Classical event-driven molecular dynamics: exact quadratic times of impact,
a binary-heap event calendar with per-body collision counters for
invalidation, and a uniform cell grid with cell-crossing events. Positions
are lazy: body i is at ``pos[i] + vel[i]*(t - tl[i])``.

Static geometry is infinite-line walls (specular or thermal) and fixed
circular posts (specular). Each post is registered in every cell that a disc
centre can occupy while touching it, so a disc only tests the posts listed
for its own cell, when it enters the cell or changes velocity.

Everything runs inside one compiled loop. Python supplies uniform random
numbers in a buffer drawn from the run's NumPy ``Generator`` so that a seed
gives the same thermal-wall draws, in the same order, as the reference
engine. The kernel returns when it reaches the requested time, the event
budget, the end of the random buffer or the event log, or an error.

Arrays travel in five tuples:

- ``B`` bodies: pos, vel, tl, rad, mass, cnt
- ``G`` grid: cell, head, nxt, prv, fp, ip
- ``S`` static geometry: wpt, wnrm, wkind, wtemp, ppos, prad, cps, cpi
  (``cps``/``cpi``: per-cell post lists in compressed-row form)
- ``C`` calendar and scalars: ht, hk, si, sf
- ``IO`` buffers and accumulators: rand, lt, lk, lf, pimp
  (``pimp``: Kahan-summed impulse on each post, [ix, cx, iy, cy] per post)

Physics matches ``core.boundaries.resolve_wall`` and
``core.contacts.resolve_elastic`` for spinless discs.
"""
from __future__ import annotations

import math
import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover - exercised only without numba
    njit = None

# Event types in the calendar.
EV_PAIR, EV_WALL, EV_CELL, EV_POST = 0, 1, 2, 3
# Wall kinds.
WALL_SPECULAR, WALL_HOT, WALL_COLD = 0, 1, 2
# Kinds in the event log.
LOG_ELASTIC, LOG_SPECULAR, LOG_HOT, LOG_COLD, LOG_POST = 0, 1, 2, 3, 4
# Return status.
DONE, MAX_EVENTS, NEED_RANDOM, LOG_FULL, HEAP_FULL = 0, 1, 2, 3, 4
ERR_OVERLAP, ERR_CELL = -1, -2
# Float scalars (sf).
SF_TIME, SF_KE = 0, 1
SF_HEAT_HOT, SF_HEAT_COLD, SF_HEAT_OTHER = 2, 4, 6   # value, correction
SF_SUPPORT_X, SF_SUPPORT_Y = 8, 10
SF_ERR_GAP = 12
SF_SIZE = 13
# Integer scalars (si).
SI_HEAP, SI_EVENTS, SI_RAND, SI_LOG, SI_LOG_ON = 0, 1, 2, 3, 4
SI_ERR_A, SI_ERR_B, SI_ERR_KIND = 5, 6, 7
SI_SIZE = 8
# Float parameters (fp) and integer parameters (ip).
FP_GEOM, FP_VEL, FP_LOX, FP_LOY, FP_CX, FP_CY = 0, 1, 2, 3, 4, 5
IP_NX, IP_NY = 0, 1


def numba_available() -> bool:
    return njit is not None


def post_cell_lists(ppos: np.ndarray, prad: np.ndarray, reach: float,
                    fp: np.ndarray, ip: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compressed-row lists of the posts each cell must test.

    A post is listed in every cell overlapping its bounding square grown by
    ``reach`` (the largest disc radius plus a margin): any disc touching the
    post has its centre in one of those cells.
    """
    nx, ny = int(ip[IP_NX]), int(ip[IP_NY])
    lists: list[list[int]] = [[] for _ in range(nx*ny)]
    for k, ((x, y), r) in enumerate(zip(ppos, prad)):
        grow = r + reach
        x0 = max(0, int(math.floor((x - grow - fp[FP_LOX])/fp[FP_CX])))
        x1 = min(nx - 1, int(math.floor((x + grow - fp[FP_LOX])/fp[FP_CX])))
        y0 = max(0, int(math.floor((y - grow - fp[FP_LOY])/fp[FP_CY])))
        y1 = min(ny - 1, int(math.floor((y + grow - fp[FP_LOY])/fp[FP_CY])))
        for cy in range(y0, y1 + 1):
            for cx in range(x0, x1 + 1):
                lists[cy*nx + cx].append(k)
    cps = np.zeros(nx*ny + 1, np.int64)
    cps[1:] = np.cumsum([len(x) for x in lists])
    cpi = np.asarray([k for x in lists for k in x], np.int64)
    return cps, cpi


if njit is not None:
    @njit(cache=True, fastmath=False, inline="always")
    def _kadd(acc, k, x):
        y = x - acc[k+1]
        t = acc[k] + y
        acc[k+1] = (t - acc[k]) - y
        acc[k] = t

    @njit(cache=True, fastmath=False, inline="always")
    def _heap_push(C, t, kind, a, b, ca, cb):
        ht, hk, si, _ = C
        i = si[SI_HEAP]
        si[SI_HEAP] = i + 1
        while i > 0:
            parent = (i - 1) >> 1
            if ht[parent] <= t:
                break
            ht[i] = ht[parent]
            for c in range(5):
                hk[i, c] = hk[parent, c]
            i = parent
        ht[i] = t
        hk[i, 0] = kind
        hk[i, 1] = a
        hk[i, 2] = b
        hk[i, 3] = ca
        hk[i, 4] = cb

    @njit(cache=True, fastmath=False)
    def _heap_pop(C):
        """Remove the root and return its time."""
        ht, hk, si, _ = C
        root_time = ht[0]
        n = si[SI_HEAP] - 1
        si[SI_HEAP] = n
        if n == 0:
            return root_time
        t = ht[n]
        k0, k1, k2, k3, k4 = hk[n, 0], hk[n, 1], hk[n, 2], hk[n, 3], hk[n, 4]
        i = 0
        while True:
            child = 2*i + 1
            if child >= n:
                break
            if child + 1 < n and ht[child+1] < ht[child]:
                child += 1
            if ht[child] >= t:
                break
            ht[i] = ht[child]
            for c in range(5):
                hk[i, c] = hk[child, c]
            i = child
        ht[i] = t
        hk[i, 0] = k0
        hk[i, 1] = k1
        hk[i, 2] = k2
        hk[i, 3] = k3
        hk[i, 4] = k4
        return root_time

    @njit(cache=True, fastmath=False, inline="always")
    def _cell_of(x, y, fp, ip):
        cx = int(math.floor((x - fp[FP_LOX]) / fp[FP_CX]))
        cy = int(math.floor((y - fp[FP_LOY]) / fp[FP_CY]))
        nx, ny = ip[IP_NX], ip[IP_NY]
        # Bodies touching the outer walls may sit a rounding error outside.
        if cx == -1 or cx == nx:
            cx = min(max(cx, 0), nx - 1)
        if cy == -1 or cy == ny:
            cy = min(max(cy, 0), ny - 1)
        if cx < 0 or cx >= nx or cy < 0 or cy >= ny:
            return -1
        return cy*nx + cx

    @njit(cache=True, fastmath=False, inline="always")
    def _grid_insert(i, c, G):
        cell, head, nxt, prv, _, _ = G
        cell[i] = c
        prv[i] = -1
        nxt[i] = head[c]
        if head[c] >= 0:
            prv[head[c]] = i
        head[c] = i

    @njit(cache=True, fastmath=False, inline="always")
    def _grid_remove(i, G):
        cell, head, nxt, prv, _, _ = G
        c = cell[i]
        if prv[i] >= 0:
            nxt[prv[i]] = nxt[i]
        else:
            head[c] = nxt[i]
        if nxt[i] >= 0:
            prv[nxt[i]] = prv[i]
        nxt[i] = -1
        prv[i] = -1

    @njit(cache=True, fastmath=False)
    def build_grid(B, G, t):
        """Place every body in its cell at time t; return -1 or a bad index."""
        pos, vel, tl, _, _, _ = B
        head, fp, ip = G[1], G[4], G[5]
        head[:] = -1
        for i in range(len(tl)):
            dt = t - tl[i]
            c = _cell_of(pos[i, 0] + vel[i, 0]*dt, pos[i, 1] + vel[i, 1]*dt, fp, ip)
            if c < 0:
                return i
            _grid_insert(i, c, G)
        return -1

    @njit(cache=True, fastmath=False, inline="always")
    def _circle_dt(dx, dy, vx, vy, r, fp):
        """Time until |(dx,dy) + (vx,vy) t| = r; inf if none, -1 if overlap.

        (dx, dy) and (vx, vy) are the relative position and velocity. Same
        touching convention as the reference disc kernel.
        """
        c = dx*dx + dy*dy - r*r
        approach = dx*vx + dy*vy
        geom, vtol = fp[FP_GEOM], fp[FP_VEL]
        if c <= geom:
            if math.sqrt(dx*dx + dy*dy) - r < -geom:
                return -1.0
            if approach >= -vtol:
                return math.inf
            return 0.0
        speed_sq = vx*vx + vy*vy
        if speed_sq <= vtol*vtol or approach >= 0.0:
            return math.inf
        disc = approach*approach - speed_sq*c
        if disc < 0.0:
            return math.inf
        q = -approach + math.sqrt(disc)
        return c/q if q > 0.0 else (-approach - math.sqrt(disc))/speed_sq

    @njit(cache=True, fastmath=False, inline="always")
    def _wall_dt(i, w, t, B, S, fp):
        pos, vel, tl, rad, _, _ = B
        wpt, wnrm = S[0], S[1]
        dt = t - tl[i]
        x, y = pos[i, 0] + vel[i, 0]*dt, pos[i, 1] + vel[i, 1]*dt
        nx, ny = wnrm[w, 0], wnrm[w, 1]
        gap = (x - wpt[w, 0])*nx + (y - wpt[w, 1])*ny - rad[i]
        vn = vel[i, 0]*nx + vel[i, 1]*ny
        if gap < -fp[FP_GEOM]:
            return -1.0, gap
        if vn >= -fp[FP_VEL]:
            return math.inf, gap
        if gap <= fp[FP_GEOM]:
            return 0.0, gap
        return gap/(-vn), gap

    @njit(cache=True, fastmath=False, inline="always")
    def _cell_dt(i, t, B, G):
        """Earliest crossing of i's own cell boundary: (dt, direction)."""
        pos, vel, tl, _, _, _ = B
        cell, fp, ip = G[0], G[4], G[5]
        nx, ny = ip[IP_NX], ip[IP_NY]
        c = cell[i]
        cx, cy = c % nx, c // nx
        dt = t - tl[i]
        x, y = pos[i, 0] + vel[i, 0]*dt, pos[i, 1] + vel[i, 1]*dt
        best, direction = math.inf, -1
        vx, vy = vel[i, 0], vel[i, 1]
        if vx > 0.0 and cx < nx - 1:
            best, direction = max(0.0, (fp[FP_LOX] + (cx+1)*fp[FP_CX] - x)/vx), 0
        elif vx < 0.0 and cx > 0:
            best, direction = max(0.0, (fp[FP_LOX] + cx*fp[FP_CX] - x)/vx), 1
        if vy > 0.0 and cy < ny - 1:
            ty = max(0.0, (fp[FP_LOY] + (cy+1)*fp[FP_CY] - y)/vy)
            if ty < best:
                best, direction = ty, 2
        elif vy < 0.0 and cy > 0:
            ty = max(0.0, (fp[FP_LOY] + cy*fp[FP_CY] - y)/vy)
            if ty < best:
                best, direction = ty, 3
        return best, direction

    @njit(cache=True, fastmath=False, inline="always")
    def _fail(C, kind, a, b, gap):
        si, sf = C[2], C[3]
        si[SI_ERR_A], si[SI_ERR_B], si[SI_ERR_KIND] = a, b, kind
        sf[SF_ERR_GAP] = gap
        return False

    @njit(cache=True, fastmath=False, inline="always")
    def _predict(i, t, only_greater, skip_cx, skip_cy, with_walls, B, G, S, C):
        """Push i's future events. Returns False on an overlap (recorded).

        only_greater: pairs only with j > i (full rebuild). skip_cx/skip_cy:
        a previous cell whose neighbourhood was already predicted (-1: none).
        with_walls: also predict the walls (not needed after a cell crossing).
        """
        pos, vel, tl, rad, _, cnt = B
        cell, head, nxt, _, fp, ip = G
        wpt, ppos, prad, cps, cpi = S[0], S[4], S[5], S[6], S[7]
        nx, ny = ip[IP_NX], ip[IP_NY]
        dti = t - tl[i]
        xi, yi = pos[i, 0] + vel[i, 0]*dti, pos[i, 1] + vel[i, 1]*dti
        if with_walls:
            for w in range(len(wpt)):
                dt, gap = _wall_dt(i, w, t, B, S, fp)
                if dt < 0.0:
                    return _fail(C, EV_WALL, i, w, gap)
                if dt < math.inf:
                    _heap_push(C, t + dt, EV_WALL, i, w, cnt[i], 0)
        c = cell[i]
        # Posts listed for this cell: tested on entry and after every
        # velocity change, so a crossing into a listed cell is never missed.
        for slot in range(cps[c], cps[c+1]):
            k = cpi[slot]
            dx, dy = ppos[k, 0] - xi, ppos[k, 1] - yi
            dt = _circle_dt(dx, dy, -vel[i, 0], -vel[i, 1], rad[i] + prad[k], fp)
            if dt < 0.0:
                return _fail(C, EV_POST, i, k, math.sqrt(dx*dx + dy*dy) - rad[i] - prad[k])
            if dt < math.inf:
                _heap_push(C, t + dt, EV_POST, i, k, cnt[i], 0)
        cx, cy = c % nx, c // nx
        for oy in range(-1, 2):
            gy = cy + oy
            if gy < 0 or gy >= ny:
                continue
            for ox in range(-1, 2):
                gx = cx + ox
                if gx < 0 or gx >= nx:
                    continue
                if skip_cx >= 0 and abs(gx - skip_cx) <= 1 and abs(gy - skip_cy) <= 1:
                    continue
                j = head[gy*nx + gx]
                while j >= 0:
                    if j != i and (not only_greater or j > i):
                        dtj = t - tl[j]
                        dx = pos[j, 0] + vel[j, 0]*dtj - xi
                        dy = pos[j, 1] + vel[j, 1]*dtj - yi
                        dt = _circle_dt(dx, dy, vel[j, 0] - vel[i, 0], vel[j, 1] - vel[i, 1],
                                        rad[i] + rad[j], fp)
                        if dt < 0.0:
                            return _fail(C, EV_PAIR, i, j,
                                         math.sqrt(dx*dx + dy*dy) - rad[i] - rad[j])
                        if dt < math.inf:
                            a, b = (i, j) if i < j else (j, i)
                            _heap_push(C, t + dt, EV_PAIR, a, b, cnt[a], cnt[b])
                    j = nxt[j]
        dt, direction = _cell_dt(i, t, B, G)
        if direction >= 0:
            _heap_push(C, t + dt, EV_CELL, i, direction, cnt[i], 0)
        return True

    @njit(cache=True, fastmath=False)
    def rebuild(B, G, S, C):
        """Recompute the whole calendar at the current time."""
        si, sf = C[2], C[3]
        si[SI_HEAP] = 0
        t = sf[SF_TIME]
        for i in range(len(B[2])):
            if not _predict(i, t, True, -1, -1, True, B, G, S, C):
                return False
        return True

    @njit(cache=True, fastmath=False, inline="always")
    def _move_to(i, t, B):
        pos, vel, tl = B[0], B[1], B[2]
        dt = t - tl[i]
        pos[i, 0] += vel[i, 0]*dt
        pos[i, 1] += vel[i, 1]*dt
        tl[i] = t

    @njit(cache=True, fastmath=False, inline="always")
    def _log(IO, si, t, kind, a, b, other, ix, iy, heat, e0, e1, mode_in, mode_out):
        lt, lk, lf = IO[1], IO[2], IO[3]
        k = si[SI_LOG]
        lt[k] = t
        lk[k, 0], lk[k, 1], lk[k, 2], lk[k, 3] = kind, a, b, other
        lf[k, 0], lf[k, 1], lf[k, 2], lf[k, 3], lf[k, 4] = ix, iy, heat, e0, e1
        lf[k, 5], lf[k, 6] = mode_in, mode_out
        si[SI_LOG] = k + 1

    @njit(cache=True, fastmath=False, inline="always")
    def _resolve_pair(a, b, t, B, C, IO):
        pos, vel, _, _, mass, cnt = B
        si, sf = C[2], C[3]
        _move_to(a, t, B)
        _move_to(b, t, B)
        dx, dy = pos[b, 0] - pos[a, 0], pos[b, 1] - pos[a, 1]
        norm = math.sqrt(dx*dx + dy*dy)
        nxv, nyv = dx/norm, dy/norm
        g = (vel[b, 0] - vel[a, 0])*nxv + (vel[b, 1] - vel[a, 1])*nyv
        e0 = sf[SF_KE]
        ix = iy = 0.0
        if g < 0.0:
            d = 1.0/mass[a] + 1.0/mass[b]
            j = -2.0*g/d
            ke0 = (0.5*mass[a]*(vel[a, 0]**2 + vel[a, 1]**2) +
                   0.5*mass[b]*(vel[b, 0]**2 + vel[b, 1]**2))
            vel[a, 0] -= j*nxv/mass[a]
            vel[a, 1] -= j*nyv/mass[a]
            vel[b, 0] += j*nxv/mass[b]
            vel[b, 1] += j*nyv/mass[b]
            ke1 = (0.5*mass[a]*(vel[a, 0]**2 + vel[a, 1]**2) +
                   0.5*mass[b]*(vel[b, 0]**2 + vel[b, 1]**2))
            sf[SF_KE] += ke1 - ke0
            ix, iy = j*nxv, j*nyv
        cnt[a] += 1
        cnt[b] += 1
        if si[SI_LOG_ON] != 0:
            _log(IO, si, t, LOG_ELASTIC, a, b, -1, ix, iy, 0.0, e0, sf[SF_KE], 0.0, 0.0)

    @njit(cache=True, fastmath=False, inline="always")
    def _resolve_wall(a, w, t, B, S, C, IO):
        vel, mass, cnt = B[1], B[4], B[5]
        wnrm, wkind, wtemp = S[1], S[2], S[3]
        si, sf = C[2], C[3]
        rand = IO[0]
        _move_to(a, t, B)
        nxw, nyw = wnrm[w, 0], wnrm[w, 1]
        vn = vel[a, 0]*nxw + vel[a, 1]*nyw   # along the inward normal
        e0 = sf[SF_KE]
        ix = iy = heat = 0.0
        mode_in = 0.5*mass[a]*vn*vn
        mode_out = mode_in
        wk = wkind[w]
        if vn < 0.0:
            m = mass[a]
            if wk == WALL_SPECULAR:
                target = -vn
            else:
                u = max(rand[si[SI_RAND]], 2.2250738585072014e-308)
                si[SI_RAND] += 1
                target = math.sqrt(-2.0*wtemp[w]*(1.0/m)*math.log(u))
            ke0 = 0.5*m*(vel[a, 0]**2 + vel[a, 1]**2)
            dv = target - vn
            vel[a, 0] += dv*nxw
            vel[a, 1] += dv*nyw
            ke1 = 0.5*m*(vel[a, 0]**2 + vel[a, 1]**2)
            sf[SF_KE] += ke1 - ke0
            mode_out = 0.5*m*target*target
            # Logged vector: impulse on the wall (reference convention).
            ix, iy = -m*dv*nxw, -m*dv*nyw
            _kadd(sf, SF_SUPPORT_X, -ix)
            _kadd(sf, SF_SUPPORT_Y, -iy)
            if wk != WALL_SPECULAR:
                heat = ke1 - ke0
                if wk == WALL_HOT:
                    _kadd(sf, SF_HEAT_HOT, heat)
                else:
                    _kadd(sf, SF_HEAT_COLD, heat)
        cnt[a] += 1
        if si[SI_LOG_ON] != 0:
            lkind = (LOG_SPECULAR if wk == WALL_SPECULAR else
                     LOG_HOT if wk == WALL_HOT else LOG_COLD)
            _log(IO, si, t, lkind, a, -1, w, ix, iy, heat, e0, sf[SF_KE], mode_in, mode_out)

    @njit(cache=True, fastmath=False, inline="always")
    def _resolve_post(a, k, t, B, S, C, IO):
        pos, vel, mass, cnt = B[0], B[1], B[4], B[5]
        ppos = S[4]
        si, sf = C[2], C[3]
        pimp = IO[4]
        _move_to(a, t, B)
        dx, dy = ppos[k, 0] - pos[a, 0], ppos[k, 1] - pos[a, 1]
        norm = math.sqrt(dx*dx + dy*dy)
        nxv, nyv = dx/norm, dy/norm          # from the disc toward the post
        vn = vel[a, 0]*nxv + vel[a, 1]*nyv
        e0 = sf[SF_KE]
        ix = iy = 0.0
        if vn > 0.0:
            m = mass[a]
            ke0 = 0.5*m*(vel[a, 0]**2 + vel[a, 1]**2)
            vel[a, 0] -= 2.0*vn*nxv
            vel[a, 1] -= 2.0*vn*nyv
            ke1 = 0.5*m*(vel[a, 0]**2 + vel[a, 1]**2)
            sf[SF_KE] += ke1 - ke0
            ix, iy = 2.0*m*vn*nxv, 2.0*m*vn*nyv   # impulse on the post
            _kadd(pimp, 4*k, ix)
            _kadd(pimp, 4*k + 2, iy)
            _kadd(sf, SF_SUPPORT_X, -ix)
            _kadd(sf, SF_SUPPORT_Y, -iy)
        cnt[a] += 1
        if si[SI_LOG_ON] != 0:
            _log(IO, si, t, LOG_POST, a, -1, k, ix, iy, 0.0, e0, sf[SF_KE], 0.0, 0.0)

    @njit(cache=True, fastmath=False)
    def advance(t_end, max_events, B, G, S, C, IO):
        """Process events up to t_end (inclusive) or until a stop condition."""
        cnt = B[5]
        cell, prv, ip = G[0], G[3], G[5]
        wkind = S[2]
        ht, hk, si, sf = C
        rand, lt = IO[0], IO[1]
        n = len(cnt)
        cap = len(ht)
        # One event replans at most two bodies against every body, wall and post.
        margin = 2*(n + len(S[0]) + len(S[5]) + 2)
        log_cap = len(lt)
        done = 0
        while True:
            if done >= max_events:
                return MAX_EVENTS
            if si[SI_HEAP] > cap - margin:
                if not rebuild(B, G, S, C):
                    return ERR_OVERLAP
                if si[SI_HEAP] > cap - margin:
                    return HEAP_FULL
            if si[SI_HEAP] == 0 or ht[0] > t_end:
                if t_end < math.inf:
                    sf[SF_TIME] = max(sf[SF_TIME], t_end)
                return DONE
            kind = hk[0, 0]
            a = hk[0, 1]
            b = hk[0, 2]
            if kind == EV_PAIR:
                valid = hk[0, 3] == cnt[a] and hk[0, 4] == cnt[b]
            else:
                valid = hk[0, 3] == cnt[a]
            if valid and kind != EV_CELL and si[SI_LOG_ON] != 0 and si[SI_LOG] >= log_cap:
                return LOG_FULL
            if (valid and kind == EV_WALL and wkind[b] != WALL_SPECULAR and
                    si[SI_RAND] >= len(rand)):
                return NEED_RANDOM
            t = _heap_pop(C)
            if not valid:
                continue
            sf[SF_TIME] = t
            if kind == EV_CELL:
                _move_to(a, t, B)
                nx = ip[IP_NX]
                old = cell[a]
                ocx, ocy = old % nx, old // nx
                ncx, ncy = ocx, ocy
                if b == 0:
                    ncx += 1
                elif b == 1:
                    ncx -= 1
                elif b == 2:
                    ncy += 1
                else:
                    ncy -= 1
                if ncx < 0 or ncx >= nx or ncy < 0 or ncy >= ip[IP_NY]:
                    si[SI_ERR_A], si[SI_ERR_B], si[SI_ERR_KIND] = a, b, EV_CELL
                    return ERR_CELL
                _grid_remove(a, G)
                _grid_insert(a, ncy*nx + ncx, G)
                if not _predict(a, t, False, ocx, ocy, False, B, G, S, C):
                    return ERR_OVERLAP
                continue
            if kind == EV_PAIR:
                _resolve_pair(a, b, t, B, C, IO)
                si[SI_EVENTS] += 1
                done += 1
                if not _predict(a, t, False, -1, -1, True, B, G, S, C):
                    return ERR_OVERLAP
                if not _predict(b, t, False, -1, -1, True, B, G, S, C):
                    return ERR_OVERLAP
                continue
            if kind == EV_WALL:
                _resolve_wall(a, b, t, B, S, C, IO)
            else:
                _resolve_post(a, b, t, B, S, C, IO)
            si[SI_EVENTS] += 1
            done += 1
            if not _predict(a, t, False, -1, -1, True, B, G, S, C):
                return ERR_OVERLAP

    @njit(cache=True, fastmath=False)
    def positions_at(t, pos, vel, tl, out):
        for i in range(len(tl)):
            dt = t - tl[i]
            out[i, 0] = pos[i, 0] + vel[i, 0]*dt
            out[i, 1] = pos[i, 1] + vel[i, 1]*dt

    @njit(cache=True, fastmath=False)
    def min_gap(t, B, G, S):
        """Smallest pair, wall or post gap at time t: (gap, kind, a, b)."""
        pos, vel, tl, rad, _, _ = B
        cell, head, nxt, _, _, ip = G
        wpt, wnrm, ppos, prad, cps, cpi = S[0], S[1], S[4], S[5], S[6], S[7]
        nx, ny = ip[IP_NX], ip[IP_NY]
        best, kind, ba, bb = math.inf, -1, -1, -1
        for i in range(len(tl)):
            dti = t - tl[i]
            xi, yi = pos[i, 0] + vel[i, 0]*dti, pos[i, 1] + vel[i, 1]*dti
            for w in range(len(wpt)):
                gap = (xi - wpt[w, 0])*wnrm[w, 0] + (yi - wpt[w, 1])*wnrm[w, 1] - rad[i]
                if gap < best:
                    best, kind, ba, bb = gap, EV_WALL, i, w
            c = cell[i]
            for slot in range(cps[c], cps[c+1]):
                k = cpi[slot]
                gap = math.hypot(ppos[k, 0] - xi, ppos[k, 1] - yi) - rad[i] - prad[k]
                if gap < best:
                    best, kind, ba, bb = gap, EV_POST, i, k
            cx, cy = c % nx, c // nx
            for oy in range(-1, 2):
                gy = cy + oy
                if gy < 0 or gy >= ny:
                    continue
                for ox in range(-1, 2):
                    gx = cx + ox
                    if gx < 0 or gx >= nx:
                        continue
                    j = head[gy*nx + gx]
                    while j >= 0:
                        if j > i:
                            dtj = t - tl[j]
                            dx = pos[j, 0] + vel[j, 0]*dtj - xi
                            dy = pos[j, 1] + vel[j, 1]*dtj - yi
                            gap = math.sqrt(dx*dx + dy*dy) - rad[i] - rad[j]
                            if gap < best:
                                best, kind, ba, bb = gap, EV_PAIR, i, j
                        j = nxt[j]
        return best, kind, ba, bb
