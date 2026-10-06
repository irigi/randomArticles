"""Event-driven kernel for discs, rotating ring hosts, fixed walls and posts.

Classical event-driven molecular dynamics: exact quadratic times of impact
where they exist, a binary-heap event calendar with per-body collision
counters for invalidation, and a uniform cell grid (cells at least one disc
diameter) with cell-crossing events. Positions and angles are lazy: body i is
at ``pos[i] + vel[i]*(t - tl[i])`` with angle ``ang[i] + om[i]*(t - tl[i])``.

Bodies are small discs (in the grid) and ring hosts (not in the grid). A host
is a *thick arc*: every point within ``cap`` of an arc of radius ``Rmid``
whose span leaves out a mouth of half-angle ``beta`` around body angle 0.
Its outer and inner surfaces are circles about the host centre, so their
contact times are quadratics that rotation does not change; only the two
rounded mouth ends (caps) move with the angle. A disc can touch a cap only
while its centre is in the wall band, which it can enter only through the
mouth, so the caps are searched by conservative advancement over that short
window alone.

Hosts meet walls, posts and other hosts as their full outer circle (the
mouth acts as a lid that only discs can pass). Contacts between a host and
anything except a wall are perfectly rough with probability ``p_rough``: the
whole relative contact velocity reverses, normal and tangential, which keeps
energy, momentum and angular momentum. Other contacts are smooth.

Each host is registered in the grid cells around it (its outer circle grown
by the largest disc radius and a margin) and re-registers after moving by
the margin. A disc tests the posts and hosts registered for its own cell when
it enters the cell or changes velocity; a host tests the discs in its cells
when its velocity changes. Any disc touching a host or post therefore sits in
one of its cells, and no contact is missed.

Python supplies uniform random numbers in a buffer drawn from the run's
NumPy ``Generator``. The kernel returns when it reaches the requested time,
the event budget, the end of the random buffer or the event log, or an error.

Arrays travel in six tuples:

- ``B`` bodies: pos, vel, tl, rad, mass, cnt, ang, om, inr (inverse inertia,
  0 for spinless bodies); ``rad`` of a host is its outer radius
- ``G`` grid: cell, head, nxt, prv, fp, ip
- ``S`` static geometry: wpt, wnrm, wkind, wtemp, ppos, prad, cps, cpi
- ``H`` hosts: hb (body of host k), bh (host of body i or -1), hg (geometry
  rows Rin, Rout, Rmid, cap, beta), reg (cell x host flags), hregc
  (registration centres), hregb (registered cell bounds x0, x1, y0, y1)
- ``C`` calendar and scalars: ht, hk, si, sf
- ``IO`` buffers and accumulators: rand, lt, lk, lf, pimp
"""
from __future__ import annotations

import math
import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover - exercised only without numba
    njit = None

# Event types in the calendar.
EV_PAIR, EV_WALL, EV_CELL, EV_POST, EV_DH, EV_HH, EV_REG = 0, 1, 2, 3, 4, 5, 6
# Wall kinds.
WALL_SPECULAR, WALL_HOT, WALL_COLD = 0, 1, 2
# Kinds in the event log.
LOG_ELASTIC, LOG_SPECULAR, LOG_HOT, LOG_COLD, LOG_POST, LOG_HOST, LOG_HOST_HOST = range(7)
# Disc-host contact features (log column "other").
FEAT_OUTER, FEAT_INNER, FEAT_CAP_PLUS, FEAT_CAP_MINUS = 0, 1, 2, 3
# Return status.
DONE, MAX_EVENTS, NEED_RANDOM, LOG_FULL, HEAP_FULL = 0, 1, 2, 3, 4
ERR_OVERLAP, ERR_CELL, ERR_CCD = -1, -2, -3
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
FP_ROUGH, FP_HMARGIN, FP_REACH = 6, 7, 8
FP_SIZE = 9
IP_NX, IP_NY = 0, 1
# Host geometry columns.
HG_RIN, HG_ROUT, HG_RMID, HG_CAP, HG_BETA = 0, 1, 2, 3, 4
# Cap search iteration limit before reporting an unresolved contact.
CAP_ITERATIONS = 200000


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
    jit = njit(cache=True, fastmath=False)
    inline = njit(cache=True, fastmath=False, inline="always")

    # -- small helpers -------------------------------------------------------

    @inline
    def _kadd(acc, k, x):
        y = x - acc[k+1]
        t = acc[k] + y
        acc[k+1] = (t - acc[k]) - y
        acc[k] = t

    @inline
    def _wrap(phi):
        return phi - 2.0*math.pi*math.floor((phi + math.pi)/(2.0*math.pi))

    @inline
    def _cross(ax, ay, bx, by):
        return ax*by - ay*bx

    @inline
    def _lazy(i, t, B):
        pos, vel, tl, ang, om = B[0], B[1], B[2], B[6], B[7]
        dt = t - tl[i]
        return pos[i, 0] + vel[i, 0]*dt, pos[i, 1] + vel[i, 1]*dt, ang[i] + om[i]*dt

    @inline
    def _ke(i, B):
        vel, mass, om, inr = B[1], B[4], B[7], B[8]
        e = 0.5*mass[i]*(vel[i, 0]**2 + vel[i, 1]**2)
        if inr[i] > 0.0:
            e += 0.5*om[i]*om[i]/inr[i]
        return e

    @inline
    def _move_to(i, t, B):
        pos, vel, tl, ang, om = B[0], B[1], B[2], B[6], B[7]
        dt = t - tl[i]
        pos[i, 0] += vel[i, 0]*dt
        pos[i, 1] += vel[i, 1]*dt
        if om[i] != 0.0:
            ang[i] = _wrap(ang[i] + om[i]*dt)
        tl[i] = t

    # -- calendar ------------------------------------------------------------

    @inline
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

    @inline
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

    @inline
    def _fail(C, kind, a, b, gap):
        si, sf = C[2], C[3]
        si[SI_ERR_A], si[SI_ERR_B], si[SI_ERR_KIND] = a, b, kind
        sf[SF_ERR_GAP] = gap
        return False

    # -- grid ----------------------------------------------------------------

    @inline
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

    @inline
    def _grid_insert(i, c, G):
        cell, head, nxt, prv, _, _ = G
        cell[i] = c
        prv[i] = -1
        nxt[i] = head[c]
        if head[c] >= 0:
            prv[head[c]] = i
        head[c] = i

    @inline
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

    @jit
    def build_grid(B, G, H, t):
        """Place every disc in its cell at time t; return -1 or a bad index."""
        cell, head, fp, ip = G[0], G[1], G[4], G[5]
        bh = H[1]
        head[:] = -1
        for i in range(len(B[2])):
            if bh[i] >= 0:
                cell[i] = -1
                continue
            x, y, _ = _lazy(i, t, B)
            c = _cell_of(x, y, fp, ip)
            if c < 0:
                return i
            _grid_insert(i, c, G)
        return -1

    @inline
    def _cell_bounds(x, y, half, fp, ip):
        nx, ny = ip[IP_NX], ip[IP_NY]
        x0 = max(0, int(math.floor((x - half - fp[FP_LOX])/fp[FP_CX])))
        x1 = min(nx - 1, int(math.floor((x + half - fp[FP_LOX])/fp[FP_CX])))
        y0 = max(0, int(math.floor((y - half - fp[FP_LOY])/fp[FP_CY])))
        y1 = min(ny - 1, int(math.floor((y + half - fp[FP_LOY])/fp[FP_CY])))
        return x0, x1, y0, y1

    @jit
    def _register(k, t, B, G, H):
        """(Re)register host k around its current centre; return old bounds."""
        fp, ip = G[4], G[5]
        hb, hg, reg, hregc, hregb = H[0], H[2], H[3], H[4], H[5]
        nx = ip[IP_NX]
        ox0, ox1, oy0, oy1 = hregb[k, 0], hregb[k, 1], hregb[k, 2], hregb[k, 3]
        if ox0 >= 0:
            for cy in range(oy0, oy1 + 1):
                for cx in range(ox0, ox1 + 1):
                    reg[cy*nx + cx, k] = 0
        x, y, _ = _lazy(hb[k], t, B)
        half = hg[k, HG_ROUT] + fp[FP_REACH] + fp[FP_HMARGIN]
        x0, x1, y0, y1 = _cell_bounds(x, y, half, fp, ip)
        for cy in range(y0, y1 + 1):
            for cx in range(x0, x1 + 1):
                reg[cy*nx + cx, k] = 1
        hregc[k, 0], hregc[k, 1] = x, y
        hregb[k, 0], hregb[k, 1], hregb[k, 2], hregb[k, 3] = x0, x1, y0, y1
        return ox0, ox1, oy0, oy1

    @jit
    def register_all(B, G, H, t):
        hregb = H[5]
        hregb[:, :] = -1
        H[3][:, :] = 0
        for k in range(len(H[0])):
            _register(k, t, B, G, H)

    # -- times of impact -----------------------------------------------------

    @inline
    def _circle_dt(dx, dy, vx, vy, r, fp):
        """Time until |(dx,dy) + (vx,vy) t| = r from outside; inf if none,
        -1 if overlapping. Same touching convention as the reference kernel."""
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

    @inline
    def _roots(dx, dy, wx, wy, level):
        """Both times with |(dx,dy) + (wx,wy) t| = level (nan, nan if none)."""
        a = wx*wx + wy*wy
        if a <= 0.0:
            return math.nan, math.nan
        b = dx*wx + dy*wy
        c = dx*dx + dy*dy - level*level
        disc = b*b - a*c
        if disc < 0.0:
            return math.nan, math.nan
        s = math.sqrt(disc)
        q = -(b + math.copysign(s, b))
        if q == 0.0:
            return 0.0, 0.0
        t1, t2 = q/a, c/q
        return (t1, t2) if t1 <= t2 else (t2, t1)

    @inline
    def _arc_point(qx, qy, theta, rmid, beta):
        """Nearest point of the host's mid-arc to q (host-centred coordinates).

        Returns (distance, feature, px, py). The arc spans every angle more
        than beta away from the mouth axis at body angle theta.
        """
        rho = math.sqrt(qx*qx + qy*qy)
        phi = _wrap(math.atan2(qy, qx) - theta)
        if abs(phi) >= beta:
            if rho <= 0.0:
                return rmid, FEAT_OUTER, rmid*math.cos(theta + beta), rmid*math.sin(theta + beta)
            return abs(rho - rmid), FEAT_OUTER, rmid*qx/rho, rmid*qy/rho
        px1, py1 = rmid*math.cos(theta + beta), rmid*math.sin(theta + beta)
        px2, py2 = rmid*math.cos(theta - beta), rmid*math.sin(theta - beta)
        d1 = math.hypot(qx - px1, qy - py1)
        d2 = math.hypot(qx - px2, qy - py2)
        if d1 <= d2:
            return d1, FEAT_CAP_PLUS, px1, py1
        return d2, FEAT_CAP_MINUS, px2, py2

    @inline
    def _in_span(dx, dy, wx, wy, th, om, tau, beta):
        qx, qy = dx + wx*tau, dy + wy*tau
        return abs(_wrap(math.atan2(qy, qx) - (th + om*tau))) >= beta

    @jit
    def _cap_search(t0, t1, dx, dy, wx, wy, th, om, rmid, beta, reach, fp):
        """First contact of a point (relative motion d + w t) with either cap.

        Conservative advancement on the distance to the cap centres, which
        move on the mid circle at angular speed om. Returns the time, inf, or
        -2 when the iteration limit is reached.
        """
        geom, vtol = fp[FP_GEOM], fp[FP_VEL]
        bound = math.sqrt(wx*wx + wy*wy) + abs(om)*rmid
        if bound <= 0.0:
            bound = 1e-300
        tau = t0
        for _ in range(CAP_ITERATIONS):
            if tau > t1:
                return math.inf
            qx, qy = dx + wx*tau, dy + wy*tau
            theta = th + om*tau
            best, cx, cy = math.inf, 0.0, 0.0
            for s in (1.0, -1.0):
                px, py = rmid*math.cos(theta + s*beta), rmid*math.sin(theta + s*beta)
                dist = math.hypot(qx - px, qy - py)
                if dist < best:
                    best, cx, cy = dist, px, py
            f = best - reach
            if f <= geom:
                # Relative velocity of the point with respect to the cap centre.
                rvx, rvy = wx + om*cy, wy - om*cx
                if (qx - cx)*rvx + (qy - cy)*rvy < -vtol*max(best, 1e-300):
                    return tau
                tau += (2.0*geom - f)/bound      # touching and separating
            else:
                tau += f/bound
        return -2.0

    @jit
    def _dh_dt(i, h, k, t, B, H, fp):
        """Time from t to the next contact of disc i with host body h (host k).

        inf if none, -1 if overlapping now, -2 if a cap search was unresolved.
        """
        vel, rad, om = B[1], B[3], B[7]
        hg = H[2]
        rin, rout, rmid, cap, beta = (hg[k, HG_RIN], hg[k, HG_ROUT], hg[k, HG_RMID],
                                      hg[k, HG_CAP], hg[k, HG_BETA])
        geom = fp[FP_GEOM]
        xi, yi, _ = _lazy(i, t, B)
        xh, yh, th = _lazy(h, t, B)
        dx, dy = xi - xh, yi - yh
        wx, wy = vel[i, 0] - vel[h, 0], vel[i, 1] - vel[h, 1]
        w = om[h]
        r = rad[i]
        dist, _, px, py = _arc_point(dx, dy, th, rmid, beta)
        if dist - cap - r < -geom:
            return -1.0
        if dist - cap - r <= geom:
            # Touching now: a contact if the disc approaches the arc point,
            # which moves with the host's rotation.
            rvx, rvy = wx + w*py, wy - w*px
            if (dx - px)*rvx + (dy - py)*rvy < -fp[FP_VEL]*dist:
                return 0.0
        lo, li = rout + r, rin - r
        rho2 = dx*dx + dy*dy
        region = 0 if rho2 >= lo*lo else (2 if rho2 <= li*li else 1)
        cursor = 0.0
        for _ in range(4):
            if region == 0:                      # outside the outer surface
                t1, t2 = _roots(dx, dy, wx, wy, lo)
                if not (t2 > cursor):            # nan or receding
                    return math.inf
                tau = max(t1, cursor)
                if _in_span(dx, dy, wx, wy, th, w, tau, beta):
                    return tau                   # outer arc
                region, cursor = 1, tau          # enters the band through the mouth
            elif region == 1:                    # in the wall band (mouth sector)
                t1i, _ = _roots(dx, dy, wx, wy, li)
                _, t2o = _roots(dx, dy, wx, wy, lo)
                inward = t1i > cursor            # false for nan
                if inward:
                    end = t1i
                elif t2o > cursor:
                    end = t2o
                elif w != 0.0:
                    end = cursor + 2.0*math.pi/abs(w)   # at rest in the band
                else:
                    return math.inf
                tau = _cap_search(cursor, end, dx, dy, wx, wy, th, w, rmid, beta,
                                  cap + r, fp)
                if tau == -2.0 or tau < math.inf:
                    return tau
                if not inward:
                    return math.inf
                region, cursor = 2, end
            else:                                # inside the cavity
                _, t2 = _roots(dx, dy, wx, wy, li)
                if not (t2 >= cursor):
                    return math.inf
                tau = t2
                if _in_span(dx, dy, wx, wy, th, w, tau, beta):
                    return tau                   # inner arc
                region, cursor = 1, tau          # leaves through the mouth
        return math.inf

    @inline
    def _wall_dt(i, w, t, B, S, fp):
        vel, rad = B[1], B[3]
        wpt, wnrm = S[0], S[1]
        x, y, _ = _lazy(i, t, B)
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

    @inline
    def _cell_dt(i, t, B, G):
        """Earliest crossing of disc i's own cell boundary: (dt, direction)."""
        vel = B[1]
        cell, fp, ip = G[0], G[4], G[5]
        nx, ny = ip[IP_NX], ip[IP_NY]
        c = cell[i]
        cx, cy = c % nx, c // nx
        x, y, _ = _lazy(i, t, B)
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

    # -- prediction ----------------------------------------------------------

    @jit
    def _push_dh(i, h, k, t, B, H, C, fp):
        cnt = B[5]
        dt = _dh_dt(i, h, k, t, B, H, fp)
        if dt == -1.0:
            return _fail(C, EV_DH, i, h, -1.0)
        if dt == -2.0:
            return _fail(C, EV_DH, i, h, math.nan)
        if dt < math.inf:
            _heap_push(C, t + dt, EV_DH, i, h, cnt[i], cnt[h])
        return True

    @jit
    def _predict_disc_hosts(i, c, t, B, H, C, fp):
        """Disc i against the hosts registered for its cell c (kept out of line)."""
        hb, reg = H[0], H[3]
        for k in range(len(hb)):
            if reg[c, k] != 0:
                if not _push_dh(i, hb[k], k, t, B, H, C, fp):
                    return False
        return True

    @inline
    def _predict_disc(i, t, only_greater, skip_cx, skip_cy, with_walls, B, G, S, H, C):
        """Push disc i's future events. Returns False on a failure (recorded).

        only_greater: disc pairs only with j > i (full rebuild). skip_cx/skip_cy:
        a previous cell whose neighbourhood was already predicted (-1: none).
        with_walls: also predict the walls (not needed after a cell crossing).
        """
        pos, vel, tl, rad, _, cnt = B[0], B[1], B[2], B[3], B[4], B[5]
        cell, head, nxt, _, fp, ip = G
        wpt, ppos, prad, cps, cpi = S[0], S[4], S[5], S[6], S[7]
        hb = H[0]
        nx, ny = ip[IP_NX], ip[IP_NY]
        xi, yi, _ = _lazy(i, t, B)
        if with_walls:
            for w in range(len(wpt)):
                dt, gap = _wall_dt(i, w, t, B, S, fp)
                if dt < 0.0:
                    return _fail(C, EV_WALL, i, w, gap)
                if dt < math.inf:
                    _heap_push(C, t + dt, EV_WALL, i, w, cnt[i], 0)
        c = cell[i]
        for slot in range(cps[c], cps[c+1]):
            k = cpi[slot]
            dx, dy = ppos[k, 0] - xi, ppos[k, 1] - yi
            dt = _circle_dt(dx, dy, -vel[i, 0], -vel[i, 1], rad[i] + prad[k], fp)
            if dt < 0.0:
                return _fail(C, EV_POST, i, k, math.sqrt(dx*dx + dy*dy) - rad[i] - prad[k])
            if dt < math.inf:
                _heap_push(C, t + dt, EV_POST, i, k, cnt[i], 0)
        if len(hb) > 0 and not _predict_disc_hosts(i, c, t, B, H, C, fp):
            return False
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

    @jit
    def _predict_host_discs(h, k, t, x0, x1, y0, y1, ox0, ox1, oy0, oy1, B, G, H, C):
        """Disc contacts for host h in cells x0..x1, y0..y1 outside the old box."""
        head, nxt, fp, ip = G[1], G[2], G[4], G[5]
        nx = ip[IP_NX]
        for cy in range(y0, y1 + 1):
            for cx in range(x0, x1 + 1):
                if ox0 <= cx <= ox1 and oy0 <= cy <= oy1:
                    continue
                j = head[cy*nx + cx]
                while j >= 0:
                    if not _push_dh(j, h, k, t, B, H, C, fp):
                        return False
                    j = nxt[j]
        return True

    @inline
    def _push_reg(h, k, t, B, G, H, C):
        vel, cnt = B[1], B[5]
        hregc = H[4]
        x, y, _ = _lazy(h, t, B)
        _, t2 = _roots(x - hregc[k, 0], y - hregc[k, 1], vel[h, 0], vel[h, 1],
                       G[4][FP_HMARGIN])
        if t2 >= 0.0:
            _heap_push(C, t + t2, EV_REG, h, k, cnt[h], 0)

    @jit
    def _predict_host(h, t, rebuild_mode, B, G, S, H, C):
        """Push host body h's events: walls, posts, hosts, discs, re-registration.

        In a rebuild, discs predict their own host contacts and host pairs are
        taken once (other host body > h).
        """
        vel, rad, cnt = B[1], B[3], B[5]
        fp = G[4]
        wpt, ppos, prad = S[0], S[4], S[5]
        hb, bh, hregb = H[0], H[1], H[5]
        k = bh[h]
        x, y, _ = _lazy(h, t, B)
        for w in range(len(wpt)):
            dt, gap = _wall_dt(h, w, t, B, S, fp)
            if dt < 0.0:
                return _fail(C, EV_WALL, h, w, gap)
            if dt < math.inf:
                _heap_push(C, t + dt, EV_WALL, h, w, cnt[h], 0)
        for p in range(len(prad)):
            dx, dy = ppos[p, 0] - x, ppos[p, 1] - y
            dt = _circle_dt(dx, dy, -vel[h, 0], -vel[h, 1], rad[h] + prad[p], fp)
            if dt < 0.0:
                return _fail(C, EV_POST, h, p, math.sqrt(dx*dx + dy*dy) - rad[h] - prad[p])
            if dt < math.inf:
                _heap_push(C, t + dt, EV_POST, h, p, cnt[h], 0)
        for k2 in range(len(hb)):
            j = hb[k2]
            if j == h or (rebuild_mode and j < h):
                continue
            xj, yj, _ = _lazy(j, t, B)
            dx, dy = xj - x, yj - y
            dt = _circle_dt(dx, dy, vel[j, 0] - vel[h, 0], vel[j, 1] - vel[h, 1],
                            rad[h] + rad[j], fp)
            if dt < 0.0:
                return _fail(C, EV_HH, h, j, math.sqrt(dx*dx + dy*dy) - rad[h] - rad[j])
            if dt < math.inf:
                a, b = (h, j) if h < j else (j, h)
                _heap_push(C, t + dt, EV_HH, a, b, cnt[a], cnt[b])
        if not rebuild_mode:
            if not _predict_host_discs(h, k, t, hregb[k, 0], hregb[k, 1], hregb[k, 2],
                                       hregb[k, 3], 1, 0, 1, 0, B, G, H, C):
                return False
        _push_reg(h, k, t, B, G, H, C)
        return True

    @inline
    def _predict(i, t, B, G, S, H, C):
        if H[1][i] < 0:
            return _predict_disc(i, t, False, -1, -1, True, B, G, S, H, C)
        return _predict_host(i, t, False, B, G, S, H, C)

    @jit
    def rebuild(B, G, S, H, C):
        """Recompute the whole calendar at the current time."""
        si, sf = C[2], C[3]
        bh = H[1]
        si[SI_HEAP] = 0
        t = sf[SF_TIME]
        for i in range(len(B[2])):
            if bh[i] < 0:
                ok = _predict_disc(i, t, True, -1, -1, True, B, G, S, H, C)
            else:
                ok = _predict_host(i, t, True, B, G, S, H, C)
            if not ok:
                return False
        return True

    # -- resolution ----------------------------------------------------------

    @inline
    def _log(IO, si, t, kind, a, b, other, ix, iy, heat, e0, e1, m5, m6):
        lt, lk, lf = IO[1], IO[2], IO[3]
        k = si[SI_LOG]
        lt[k] = t
        lk[k, 0], lk[k, 1], lk[k, 2], lk[k, 3] = kind, a, b, other
        lf[k, 0], lf[k, 1], lf[k, 2], lf[k, 3], lf[k, 4] = ix, iy, heat, e0, e1
        lf[k, 5], lf[k, 6] = m5, m6
        si[SI_LOG] = k + 1

    @inline
    def _rough_draw(C, IO, fp):
        p = fp[FP_ROUGH]
        if p >= 1.0:
            return True
        if p <= 0.0:
            return False
        si, rand = C[2], IO[0]
        u = rand[si[SI_RAND]]
        si[SI_RAND] += 1
        return u < p

    @jit
    def _contact_impulse(a, b, cx, cy, nx, ny, rough, B):
        """Elastic impulse between bodies a and b (b = -1: fixed support).

        The normal points from a toward b. Rough: the whole relative contact
        velocity reverses; smooth: only its normal part. The impulse J acts
        on b and -J on a. Returns (Jx, Jy, applied).

        Body a is always a circle touched along its own normal, so its lever
        arm is exactly rad[a]*n: no rounding torque in smooth contacts.
        """
        pos, vel, rad, mass, om, inr = B[0], B[1], B[3], B[4], B[7], B[8]
        tx, ty = -ny, nx
        rax, ray = rad[a]*nx, rad[a]*ny
        uax = vel[a, 0] - om[a]*ray
        uay = vel[a, 1] + om[a]*rax
        ima, ia = 1.0/mass[a], inr[a]
        if b >= 0:
            rbx, rby = cx - pos[b, 0], cy - pos[b, 1]
            ubx = vel[b, 0] - om[b]*rby
            uby = vel[b, 1] + om[b]*rbx
            imb, ib = 1.0/mass[b], inr[b]
        else:
            rbx = rby = ubx = uby = imb = ib = 0.0
        ux, uy = ubx - uax, uby - uay
        gn, gt = ux*nx + uy*ny, ux*tx + uy*ty
        if gn >= 0.0:
            return 0.0, 0.0, False
        an, at = 0.0, rad[a]
        bn, bt = _cross(rbx, rby, nx, ny), _cross(rbx, rby, tx, ty)
        dnn = ima + imb + ia*an*an + ib*bn*bn
        if rough:
            dtt = ima + imb + ia*at*at + ib*bt*bt
            dnt = ia*an*at + ib*bn*bt
            det = dnn*dtt - dnt*dnt
            jn = -2.0*(dtt*gn - dnt*gt)/det
            jt = -2.0*(dnn*gt - dnt*gn)/det
        else:
            jn, jt = -2.0*gn/dnn, 0.0
        jx, jy = jn*nx + jt*tx, jn*ny + jt*ty
        vel[a, 0] -= jx*ima
        vel[a, 1] -= jy*ima
        om[a] -= ia*(an*jn + at*jt)          # r x J in the (n, t) basis
        if b >= 0:
            vel[b, 0] += jx*imb
            vel[b, 1] += jy*imb
            om[b] += ib*(bn*jn + bt*jt)
        return jx, jy, True

    @inline
    def _resolve_pair(a, b, t, B, C, IO):
        pos, vel, mass, cnt = B[0], B[1], B[4], B[5]
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
            ke0 = _ke(a, B) + _ke(b, B)
            vel[a, 0] -= j*nxv/mass[a]
            vel[a, 1] -= j*nyv/mass[a]
            vel[b, 0] += j*nxv/mass[b]
            vel[b, 1] += j*nyv/mass[b]
            sf[SF_KE] += _ke(a, B) + _ke(b, B) - ke0
            ix, iy = j*nxv, j*nyv
        cnt[a] += 1
        cnt[b] += 1
        if si[SI_LOG_ON] != 0:
            _log(IO, si, t, LOG_ELASTIC, a, b, -1, ix, iy, 0.0, e0, sf[SF_KE], 0.0, 0.0)

    @inline
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

    @inline
    def _resolve_post(a, k, t, rough, B, S, C, IO):
        """Body a meets post k: smooth for discs, rough with probability for hosts."""
        pos, rad, cnt = B[0], B[3], B[5]
        ppos = S[4]
        si, sf = C[2], C[3]
        pimp = IO[4]
        _move_to(a, t, B)
        dx, dy = ppos[k, 0] - pos[a, 0], ppos[k, 1] - pos[a, 1]
        norm = math.sqrt(dx*dx + dy*dy)
        nxv, nyv = dx/norm, dy/norm          # from the body toward the post
        e0 = sf[SF_KE]
        ke0 = _ke(a, B)
        # J acts on the post (b = -1), -J on the body.
        ix, iy, applied = _contact_impulse(a, -1, pos[a, 0] + rad[a]*nxv,
                                           pos[a, 1] + rad[a]*nyv, nxv, nyv, rough, B)
        if applied:
            sf[SF_KE] += _ke(a, B) - ke0
            _kadd(pimp, 4*k, ix)
            _kadd(pimp, 4*k + 2, iy)
            _kadd(sf, SF_SUPPORT_X, -ix)
            _kadd(sf, SF_SUPPORT_Y, -iy)
        cnt[a] += 1
        if si[SI_LOG_ON] != 0:
            _log(IO, si, t, LOG_POST, a, -1, k, ix, iy, 0.0, e0, sf[SF_KE],
                 1.0 if rough else 0.0, 0.0)

    @jit
    def _resolve_dh(i, h, t, rough, B, H, C, IO):
        """Disc i meets host body h at the nearest point of the thick arc."""
        pos, rad, cnt = B[0], B[3], B[5]
        bh, hg = H[1], H[2]
        si, sf = C[2], C[3]
        k = bh[h]
        _move_to(i, t, B)
        _move_to(h, t, B)
        qx, qy = pos[i, 0] - pos[h, 0], pos[i, 1] - pos[h, 1]
        dist, feat, px, py = _arc_point(qx, qy, B[6][h], hg[k, HG_RMID], hg[k, HG_BETA])
        if feat == FEAT_OUTER and qx*qx + qy*qy < hg[k, HG_RMID]**2:
            feat = FEAT_INNER
        nxv, nyv = (px - qx)/dist, (py - qy)/dist     # from the disc toward the arc
        e0 = sf[SF_KE]
        ke0 = _ke(i, B) + _ke(h, B)
        jx, jy, applied = _contact_impulse(i, h, pos[i, 0] + rad[i]*nxv,
                                           pos[i, 1] + rad[i]*nyv, nxv, nyv, rough, B)
        if applied:
            sf[SF_KE] += _ke(i, B) + _ke(h, B) - ke0
        cnt[i] += 1
        cnt[h] += 1
        if si[SI_LOG_ON] != 0:
            _log(IO, si, t, LOG_HOST, i, h, feat, jx, jy, 0.0, e0, sf[SF_KE],
                 1.0 if rough else 0.0, 0.0)

    @jit
    def _resolve_hh(a, b, t, rough, B, C, IO):
        """Hosts meet as their outer circles."""
        pos, rad, cnt = B[0], B[3], B[5]
        si, sf = C[2], C[3]
        _move_to(a, t, B)
        _move_to(b, t, B)
        dx, dy = pos[b, 0] - pos[a, 0], pos[b, 1] - pos[a, 1]
        norm = math.sqrt(dx*dx + dy*dy)
        nxv, nyv = dx/norm, dy/norm
        e0 = sf[SF_KE]
        ke0 = _ke(a, B) + _ke(b, B)
        jx, jy, applied = _contact_impulse(a, b, pos[a, 0] + rad[a]*nxv,
                                           pos[a, 1] + rad[a]*nyv, nxv, nyv, rough, B)
        if applied:
            sf[SF_KE] += _ke(a, B) + _ke(b, B) - ke0
        cnt[a] += 1
        cnt[b] += 1
        if si[SI_LOG_ON] != 0:
            _log(IO, si, t, LOG_HOST_HOST, a, b, -1, jx, jy, 0.0, e0, sf[SF_KE],
                 1.0 if rough else 0.0, 0.0)

    # -- main loop -----------------------------------------------------------

    @jit
    def advance(t_end, max_events, B, G, S, H, C, IO):
        """Process events up to t_end (inclusive) or until a stop condition."""
        cnt = B[5]
        cell, ip, fp = G[0], G[5], G[4]
        wkind = S[2]
        bh, hregb = H[1], H[5]
        ht, hk, si, sf = C
        rand, lt = IO[0], IO[1]
        n = len(cnt)
        cap = len(ht)
        # One event replans at most two bodies against every body, wall and post.
        margin = 2*(n + len(S[0]) + len(S[5]) + len(H[0]) + 4)
        log_cap = len(lt)
        draws_rough = 0.0 < fp[FP_ROUGH] < 1.0
        done = 0
        while True:
            if done >= max_events:
                return MAX_EVENTS
            if si[SI_HEAP] > cap - margin:
                if not rebuild(B, G, S, H, C):
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
            if kind == EV_PAIR or kind == EV_DH or kind == EV_HH:
                valid = hk[0, 3] == cnt[a] and hk[0, 4] == cnt[b]
            else:
                valid = hk[0, 3] == cnt[a]
            if valid:
                if (kind != EV_CELL and kind != EV_REG and si[SI_LOG_ON] != 0 and
                        si[SI_LOG] >= log_cap):
                    return LOG_FULL
                needs = ((kind == EV_WALL and wkind[b] != WALL_SPECULAR) or
                         (draws_rough and (kind == EV_DH or kind == EV_HH or
                                           (kind == EV_POST and bh[a] >= 0))))
                if needs and si[SI_RAND] >= len(rand):
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
                if not _predict_disc(a, t, False, ocx, ocy, False, B, G, S, H, C):
                    return ERR_CCD if math.isnan(sf[SF_ERR_GAP]) else ERR_OVERLAP
                continue
            if kind == EV_REG:
                k = b
                ox0, ox1, oy0, oy1 = _register(k, t, B, G, H)
                if not _predict_host_discs(a, k, t, hregb[k, 0], hregb[k, 1], hregb[k, 2],
                                           hregb[k, 3], ox0, ox1, oy0, oy1, B, G, H, C):
                    return ERR_CCD if math.isnan(sf[SF_ERR_GAP]) else ERR_OVERLAP
                _push_reg(a, k, t, B, G, H, C)
                continue
            if kind == EV_PAIR:
                _resolve_pair(a, b, t, B, C, IO)
            elif kind == EV_WALL:
                _resolve_wall(a, b, t, B, S, C, IO)
            elif kind == EV_POST:
                _resolve_post(a, b, t, bh[a] >= 0 and _rough_draw(C, IO, fp), B, S, C, IO)
            elif kind == EV_DH:
                _resolve_dh(a, b, t, _rough_draw(C, IO, fp), B, H, C, IO)
            else:
                _resolve_hh(a, b, t, _rough_draw(C, IO, fp), B, C, IO)
            si[SI_EVENTS] += 1
            done += 1
            if not _predict(a, t, B, G, S, H, C):
                return ERR_CCD if math.isnan(sf[SF_ERR_GAP]) else ERR_OVERLAP
            if kind == EV_PAIR or kind == EV_DH or kind == EV_HH:
                if not _predict(b, t, B, G, S, H, C):
                    return ERR_CCD if math.isnan(sf[SF_ERR_GAP]) else ERR_OVERLAP

    # -- sampling ------------------------------------------------------------

    @jit
    def state_at(t, B, out_pos, out_ang):
        pos, vel, tl, ang, om = B[0], B[1], B[2], B[6], B[7]
        for i in range(len(tl)):
            dt = t - tl[i]
            out_pos[i, 0] = pos[i, 0] + vel[i, 0]*dt
            out_pos[i, 1] = pos[i, 1] + vel[i, 1]*dt
            out_ang[i] = _wrap(ang[i] + om[i]*dt)

    @jit
    def min_gap(t, B, G, S, H):
        """Smallest gap at time t: (gap, event kind, a, b)."""
        rad = B[3]
        cell, head, nxt, _, _, ip = G
        wpt, wnrm, ppos, prad, cps, cpi = S[0], S[1], S[4], S[5], S[6], S[7]
        hb, bh, hg, reg = H[0], H[1], H[2], H[3]
        nx, ny = ip[IP_NX], ip[IP_NY]
        best, kind, ba, bb = math.inf, -1, -1, -1
        for i in range(len(rad)):
            xi, yi, _ = _lazy(i, t, B)
            for w in range(len(wpt)):
                gap = (xi - wpt[w, 0])*wnrm[w, 0] + (yi - wpt[w, 1])*wnrm[w, 1] - rad[i]
                if gap < best:
                    best, kind, ba, bb = gap, EV_WALL, i, w
            if bh[i] >= 0:
                for p in range(len(prad)):
                    gap = math.hypot(ppos[p, 0] - xi, ppos[p, 1] - yi) - rad[i] - prad[p]
                    if gap < best:
                        best, kind, ba, bb = gap, EV_POST, i, p
                for k2 in range(len(hb)):
                    j = hb[k2]
                    if j > i:
                        xj, yj, _ = _lazy(j, t, B)
                        gap = math.hypot(xj - xi, yj - yi) - rad[i] - rad[j]
                        if gap < best:
                            best, kind, ba, bb = gap, EV_HH, i, j
                continue
            c = cell[i]
            for slot in range(cps[c], cps[c+1]):
                p = cpi[slot]
                gap = math.hypot(ppos[p, 0] - xi, ppos[p, 1] - yi) - rad[i] - prad[p]
                if gap < best:
                    best, kind, ba, bb = gap, EV_POST, i, p
            for k in range(len(hb)):
                if reg[c, k] != 0:
                    h = hb[k]
                    xh, yh, th = _lazy(h, t, B)
                    dist, _, _, _ = _arc_point(xi - xh, yi - yh, th, hg[k, HG_RMID],
                                               hg[k, HG_BETA])
                    gap = dist - hg[k, HG_CAP] - rad[i]
                    if gap < best:
                        best, kind, ba, bb = gap, EV_DH, i, h
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
                            xj, yj, _ = _lazy(j, t, B)
                            gap = math.hypot(xj - xi, yj - yi) - rad[i] - rad[j]
                            if gap < best:
                                best, kind, ba, bb = gap, EV_PAIR, i, j
                        j = nxt[j]
        return best, kind, ba, bb
