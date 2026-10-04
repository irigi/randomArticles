"""Session 3: independently reproduce value of perfect switch-time information.

REQUIRES the checkpointed multistart known-time N=480 solutions for integer
S=22..60 and half-integer S=26.5,27.5,28.5,29.5. Does *not* solve the NLP
implicitly. Generate/rechallenge the points independently with:
    python numerics/known_grid_checkpoint.py --start 22 --end 60 --step 1
    python numerics/known_grid_checkpoint.py --start 22 --end 60 --step 1 --reverse
    python numerics/known_grid_checkpoint.py --times 26.5 27.5 28.5 29.5

The original reproduction script wrongly used Ts>=26 without Ts<=30 when
averaging U(26,30). That bug is removed. Its claimed 0.5 s sampling also
actually used 1 s; this version records the actual nodes and quadrature.

At early S <= D0/VMAX, the known-time optimum is the unconstrained free-flow
route. Beyond that threshold, use monotone piecewise cubic interpolation of
computed local-best costs, including the exact threshold knot. This is a
numerical approximation *without* a global-optimality certificate.

Regret R(s) = J_feedback(s) - J_known(s), where the feedback realization uses
transcription-cell running-cost rates and linearly interpolated switch states.
Each prior's expectation is checked against the corresponding original
survival-weighted NLP objective via Gauss-Legendre quadrature cell by cell.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from numpy.polynomial.legendre import leggauss
from model import Model, D0, VMAX

ROOT = Path(__file__).resolve().parents[1]
SWEEP = ROOT / 'diagnostics' / 'known_sweep'
DATA = ROOT / 'data'
MODEL = Model(300, rest_mode='balance')
GRID_CELLS = 480
HALVES = [26.5, 27.5, 28.5, 29.5]


def load_known():
    base = MODEL.base_route_cost()
    knots = np.r_[np.arange(0., 22.), D0/VMAX,
                  np.arange(22., 61.), HALVES]
    cost = []
    for s in knots:
        if s <= D0/VMAX + 1e-12:
            cost.append(base)
        else:
            path = SWEEP / f'balance_{s:06.2f}_N{GRID_CELLS}.npz'
            if not path.is_file():
                raise FileNotFoundError(f'Missing verified known-time solve {path}; run known_grid_checkpoint.py')
            with np.load(path, allow_pickle=False) as f:
                cost.append(float(f['cost']))
    order = np.argsort(knots)
    knots, cost = knots[order], np.asarray(cost)[order]
    assert np.all(np.diff(knots)>0) and np.all(np.diff(cost)>-1e-5)
    return knots, cost, PchipInterpolator(knots, cost, extrapolate=False)


def load_policy(name):
    with np.load(DATA / name, allow_pickle=False) as z:
        return {k:z[k] for k in z.files}


def feedback_realized(policy, switch_s):
    """Conditional cost at arbitrary switch times within a stored red branch.

    Running cost is taken constant at each NLP cell's midpoint (the same
    quadrature as the objective); the within-cell switch state is linearly
    interpolated. This intentionally avoids the older uncapped and unvalidated
    extrapolation at the horizon endpoint.
    """
    s = np.asarray(switch_s, dtype=float)
    t, v, d, u = (policy[k] for k in ('t','v','d','u'))
    if np.any((s<t[0]-1e-10)|(s>t[-1]+1e-10)):
        raise ValueError('switch outside conditional-red horizon')
    dt = np.diff(t)
    running = np.asarray(MODEL.running_cost(.5*(v[:-1]+v[1:]),u))
    running_cumulative = np.r_[0,np.cumsum(dt*running)]
    cell = np.clip(np.searchsorted(t,s,side='right')-1,0,len(u)-1)
    alpha = np.clip((s-t[cell])/dt[cell],0,1)
    at_green_v = v[cell]*(1-alpha) + v[cell+1]*alpha
    at_green_d = d[cell]*(1-alpha) + d[cell+1]*alpha
    return running_cumulative[cell] + (s-t[cell])*running[cell] + MODEL.green_value(at_green_d,at_green_v)


def cell_gauss_integral(policy, f, lo, hi, order=8):
    """Uniform-prior expectation evaluated within each policy time cell."""
    t=policy['t']; left=np.maximum(t[:-1],lo); right=np.minimum(t[1:],hi)
    keep=right>left
    nodes,weights=leggauss(order)
    xx=(left[keep,None]+right[keep,None])/2 + (right[keep,None]-left[keep,None])/2*nodes
    ww=(right[keep,None]-left[keep,None])/2*weights
    return float(np.sum(ww*f(xx)))/(hi-lo)


def calculate(write=True):
    knots,cost,known=load_known()
    T=np.arange(0.,60.000001,.25)
    known_dense=known(T)
    results={}
    for label,filename,lo,hi in (('U_0_60','u_0_60.npz',0.,60.),
                                  ('U_26_30','u_26_30.npz',26.,30.)):
        policy=load_policy(filename)
        stored=float(policy['cost'])
        fback=lambda s: feedback_realized(policy,s)
        expected_feedback=cell_gauss_integral(policy,fback,lo,hi)
        # From the PCHIP antiderivative: actual interpolation integral, NOT
        # a mask applied to a grid extending beyond the prior support.
        perfect=float(known.integrate(lo,hi))/(hi-lo)
        independent_regret=cell_gauss_integral(policy,lambda s:fback(s)-known(s),lo,hi)
        x=np.arange(lo,hi+.001,.25)
        r=fback(x)-known(x)
        results[label]={
            'support_s':[lo,hi],
            'known_reference_cells':GRID_CELLS,
            'perfect_information_expected_CZK':perfect,
            'feedback_nlp_expected_CZK':stored,
            'feedback_realized_integrated_CZK':expected_feedback,
            'feedback_quadrature_residual_CZK':expected_feedback-stored,
            'value_information_using_nlp_CZK':stored-perfect,
            'integrated_regret_CZK':independent_regret,
            'minimum_grid_regret_CZK':float(min(r)),
            'maximum_grid_regret_CZK':float(max(r)),
        }
        assert abs(expected_feedback-stored)<2e-4, (label,'reproduction discrepancy')
        assert np.min(r)>-4e-4, (label,'unexpected negative regret')
        results[label]['sensitivity_note'] = ('Known-time reference is the best feasible local '
            'solution found at 480 transcription cells; PCHIP interpolation/integration '
            'is not a global optimum or temporal-discretization certificate.')
    if write:
        DATA.mkdir(exist_ok=True)
        broad=load_policy('u_0_60.npz');late=load_policy('u_26_30.npz')
        ltime=np.arange(26.,30.000001,.025)
        np.savez(DATA/'voi_grid.npz', T=knots, C=cost, grid_cells=GRID_CELLS,
                 broad_pk=results['U_0_60']['perfect_information_expected_CZK'],
                 late_pk=results['U_26_30']['perfect_information_expected_CZK'])
        np.savez(DATA/'voi_regret.npz', S=T, known=known_dense,
                 feedback_broad=feedback_realized(broad,T),
                 regret_broad=feedback_realized(broad,T)-known_dense,
                 late_S=ltime, late_known=known(ltime),
                 late_feedback=feedback_realized(late,ltime),
                 late_regret=feedback_realized(late,ltime)-known(ltime))
        (ROOT/'diagnostics'/'session3_voi.json').write_text(json.dumps(results,indent=2)+'\n')
    return results


if __name__=='__main__':
    out=calculate()
    for k,v in out.items():
        print(f"{k}: known={v['perfect_information_expected_CZK']:.9f}, "
              f"feedback={v['feedback_nlp_expected_CZK']:.9f}, "
              f"VOI={v['value_information_using_nlp_CZK']:.9f}, "
              f"quadrature residual={v['feedback_quadrature_residual_CZK']:.3g}")
