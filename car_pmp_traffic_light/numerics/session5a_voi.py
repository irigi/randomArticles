"""Recompute capped-mechanics VOI and regret, independent of original-model caches.

Original (uncapped) VOI calculations are never overwritten.  The feedback
realization uses the same convention and quadrature as session3 but a capped
Model for both running and green-continuation cost.
"""
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from numpy.polynomial.legendre import leggauss
from model import Model, D0, VMAX, published_accel_cap

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'
KNOWN=ROOT/'diagnostics'/'capped_known_sweep'
MODEL=Model(300,accel_cap=published_accel_cap(),brake_a=4,rest_mode='balance')
N=480
HALVES=(26.5,27.5,28.5,29.5)

def load_known():
    times=np.unique(np.r_[np.arange(0.,22.),D0/VMAX, np.arange(22.,61.),HALVES])
    cost=[]
    for s in times:
        if s <= D0/VMAX+1e-12:
            cost.append(MODEL.base_route_cost());continue
        f=KNOWN/f'capped_b4_{s:06.2f}_N{N}.npz'
        if not f.exists():raise FileNotFoundError(f)
        with np.load(f,allow_pickle=False) as z:
            assert str(z['model_id'])=='capped_b4_balance'
            assert np.isclose(z['accel_cap'],MODEL.accel_cap)
            assert np.isclose(z['brake_a'],MODEL.brake_a)
            cost.append(float(z['cost']))
    cost=np.asarray(cost)
    assert np.all(np.isfinite(cost))
    assert np.min(np.diff(cost))>-1e-4, ('nonmonotone known S',times,cost)
    return times,cost,PchipInterpolator(times,cost,extrapolate=False)

def policy(name):
    with np.load(DATA/name,allow_pickle=False) as z:
        return {k:z[k] for k in z.files}

def feedback_realized(p,s):
    s=np.asarray(s,float)
    t,v,d,u=(p[k] for k in ('t','v','d','u'))
    if np.any(s<t[0]-1e-10) or np.any(s>t[-1]+1e-10):
        raise ValueError('outside conditional-red horizon')
    dt=np.diff(t)
    running=np.asarray(MODEL.running_cost(.5*(v[:-1]+v[1:]),u))
    cumul=np.r_[0,np.cumsum(dt*running)]
    i=np.clip(np.searchsorted(t,s,side='right')-1,0,len(u)-1)
    alpha=np.clip((s-t[i])/dt[i],0,1)
    return cumul[i]+(s-t[i])*running[i]+MODEL.green_value(d[i]*(1-alpha)+d[i+1]*alpha,v[i]*(1-alpha)+v[i+1]*alpha)

def integrate_cells(p,f,lo,hi,order=8):
    t=p['t']; left=np.maximum(t[:-1],lo);right=np.minimum(t[1:],hi)
    use=right>left
    x,w=leggauss(order)
    xx=(left[use,None]+right[use,None])/2+(right[use,None]-left[use,None])/2*x
    ww=(right[use,None]-left[use,None])/2*w
    return float(np.sum(ww*f(xx)))/(hi-lo)

def calculate(write=True):
    x,y,interp=load_known()
    results={}
    for name,file,lo,hi in [('U_0_60','u_0_60_capped.npz',0.,60.),('U_26_30','u_26_30_capped.npz',26.,30.)]:
        p=policy(file)
        assert str(p['model_id'])=='capped_b4_balance'
        assert np.isclose(p['accel_cap'],MODEL.accel_cap) and np.isclose(p['brake_a'],MODEL.brake_a)
        f=lambda s: feedback_realized(p,s)
        integrated=integrate_cells(p,f,lo,hi)
        perfect=float(interp.integrate(lo,hi)/(hi-lo))
        regret=integrate_cells(p,lambda s: f(s)-interp(s),lo,hi)
        grid=np.linspace(lo,hi,1+round((hi-lo)*4))
        r=f(grid)-interp(grid)
        results[name]={'model':'capped_b4_balance','support_s':[lo,hi],
         'known_reference_cells':N,'feedback_cells':len(p['u']),
         'feedback_nlp_expected_CZK':float(p['cost']),
         'feedback_independent_integration_CZK':integrated,
         'quadrature_residual_CZK':integrated-float(p['cost']),
         'perfect_information_expected_CZK':perfect,
         'value_information_CZK':float(p['cost'])-perfect,
         'integrated_regret_CZK':regret,
         'regret_range_CZK':[float(min(r)),float(max(r))]}
        assert abs(results[name]['quadrature_residual_CZK'])<2e-4,results[name]
        assert min(r)>-5e-4,results[name]
    if write:
        (ROOT/'diagnostics'/'session5a_voi.json').write_text(json.dumps(results,indent=2)+'\n')
        np.savez_compressed(DATA/'voi_capped_grid.npz',S=x,cost=y,model_id='capped_b4_balance',cells=N)
        bb=policy('u_0_60_capped.npz');ll=policy('u_26_30_capped.npz')
        ts=np.arange(0,60.001,.25);late_ts=np.arange(26,30.001,.025)
        np.savez_compressed(DATA/'voi_capped_regret.npz',S=ts,known=interp(ts),
         feedback_broad=feedback_realized(bb,ts),regret_broad=feedback_realized(bb,ts)-interp(ts),
         late_S=late_ts,late_known=interp(late_ts),late_feedback=feedback_realized(ll,late_ts),
         late_regret=feedback_realized(ll,late_ts)-interp(late_ts),model_id='capped_b4_balance')
    return results
if __name__=='__main__':
    print(json.dumps(calculate(),indent=2))
