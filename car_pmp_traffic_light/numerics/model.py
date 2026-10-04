"""Vehicle, fuel, and green-continuation model used by the article.

The baseline model reproduces the companion article.  A sensitivity variant can
cap net full-throttle acceleration and lower the full-brake deceleration.
"""
from __future__ import annotations
import numpy as np
import casadi as ca
from scipy.integrate import cumulative_trapezoid
from numpy.polynomial import Chebyshev

# Geometry and free-flow speed
D0 = 300.0
POST = 100.0
VMAX = 50.0 / 3.6

# Fuel fit, Q in L/h with v in m/s
Q0 = 0.6
C1 = 0.120713517724
C2 = -0.000505754823921
C3 = 0.0000881390936848
FUEL_PRICE = 32.0
QI = Q0 / 3600.0

# Mechanics inherited from the companion article
R0 = 0.11772
KD = 2.712966667e-4
A_POWER = 66.0
VREG = 1.0
BASE_BRAKE = 7.0


def fuel_cruise_lph(v):
    return Q0 + C1*v + C2*v*v + C3*v*v*v


def road_load(v):
    return R0 + KD*v*v


class Model:
    def __init__(self, time_cost_czk_h=300.0, accel_cap=None,
                 brake_a=BASE_BRAKE, rest_mode="balance"):
        """rest_mode: balance (legacy moving ODE even at v=0), or
        contact (static ground reaction prevents negative velocity at rest).
        In contact mode u=b=0 holds the vehicle stationary without engine
        power. Use contact mode for physical stop/wait tests; results from
        different rest conventions should be reported separately.
        """
        if rest_mode not in ("balance", "contact"):
            raise ValueError("rest_mode must be balance or contact")
        self.rest_mode = rest_mode
        self.time_cost_czk_h = float(time_cost_czk_h)
        self.ct = self.time_cost_czk_h / 3600.0
        self.accel_cap = accel_cap
        self.brake_a = float(brake_a)
        self._build_green_value()

    def throttle_accel(self, v):
        """Acceleration coefficient multiplying throttle before road load."""
        raw = A_POWER / (v + VREG)
        if self.accel_cap is None:
            return raw
        cap_before_road = self.accel_cap + road_load(v)
        if isinstance(v, (ca.MX, ca.SX)):
            return ca.fmin(raw, cap_before_road)
        return np.minimum(raw, cap_before_road)

    def net_full_accel(self, v):
        return self.throttle_accel(v) - road_load(v)

    def cruise_throttle(self, v):
        return road_load(v) / self.throttle_accel(v)

    def gamma(self, v):
        # Chosen so that steady cruise reproduces the fitted Q(v) exactly.
        return (fuel_cruise_lph(v)/3600.0 - QI) / self.cruise_throttle(v)

    def running_cost(self, v, u):
        return self.ct + FUEL_PRICE * (QI + self.gamma(v)*u)

    def _build_green_value(self):
        self.g0 = (self.ct + FUEL_PRICE*fuel_cruise_lph(VMAX)/3600.0) / VMAX
        vv = np.linspace(0.0, VMAX, 40001)
        hp = -(self.ct + FUEL_PRICE*(QI + self.gamma(vv)) - self.g0*vv) / self.net_full_accel(vv)
        hh = cumulative_trapezoid(hp[::-1], vv[::-1], initial=0.0)[::-1]
        self._pH = Chebyshev.fit(vv[::200], hh[::200], 14, domain=(0.0, VMAX))
        self._coeff = self._pH.coef

    def H(self, v):
        z = 2*v/VMAX - 1
        if isinstance(v, (ca.MX, ca.SX)):
            t0 = 1.0
            value = self._coeff[0]
            t1 = z
            value += self._coeff[1]*t1
            for i in range(2, len(self._coeff)):
                t2 = 2*z*t1 - t0
                value += self._coeff[i]*t2
                t0, t1 = t1, t2
            return value
        return self._pH(v)

    def green_value(self, d, v):
        return self.g0*(d + POST) + self.H(v)

    def rhs_v(self, v, u, b):
        raw = self.throttle_accel(v)*u - road_load(v) - self.brake_a*b
        if self.rest_mode == 'balance':
            # Legacy smooth ODE.  At v=0 a stationary solution satisfies
            # a_T(0)*u=r(0)+a_b*b, so u=b=0 is NOT stationary.
            return raw
        # Unilateral static-ground contact at v=0: negative acceleration is
        # balanced by an implicit normal/rolling contact reaction.  A stopped
        # car with u=b=0 remains still; positive tractive force can move it.
        # This kink is a nonsmooth physics boundary, and NLP convergence is
        # assessed separately from the legacy model in the test suite.
        if isinstance(v, (ca.MX, ca.SX)):
            return ca.if_else(v <= 0, ca.fmax(raw, 0), raw)
        return np.where(np.asarray(v) <= 0, np.maximum(raw, 0), raw)

    def base_route_cost(self):
        return self.g0*(D0 + POST)


def published_accel_cap():
    """Net-acceleration cap calibrated so this simplified model gives 0-100 km/h in 8.4 s.

    Numerically precomputed from the published Octavia 2.0 TDI 110 kW time.
    """
    return 4.136672736266082
