"""Fast reproducible invariants, no IPOPT.  Run after capped sweep and VOI."""
import sys,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'numerics'))
from model import Model,published_accel_cap,D0,VMAX
from session5a_voi import load_known, policy,feedback_realized,integrate_cells,calculate

m=Model(300,accel_cap=published_accel_cap(),brake_a=4,rest_mode='balance')
assert np.isclose(m.net_full_accel(0),published_accel_cap())
x,y,f=load_known()
assert len(x)==66, len(x)  # 22 free-flow integer knots, threshold, 39 later integers, four half-second knots
assert np.isclose(x[0],0.)
assert np.all(np.diff(x)>0)
assert np.min(np.diff(y))>=-1e-4
assert np.isclose(float(f(D0/VMAX)),m.base_route_cost())
for s in (28,40,50,60):
    fn=ROOT/'diagnostics'/'capped_known_sweep'/f'capped_b4_{s:06.2f}_N480.npz'
    with np.load(fn,allow_pickle=False) as z:
        assert str(z['model_id'])=='capped_b4_balance'
        assert np.isclose(float(f(s)),float(z['cost']))
        assert z['t'].size==481
        assert np.min(z['d'])>=-2e-6 and np.min(z['v'])>=-2e-6
        assert np.max(z['v'])<=VMAX+2e-6
for name,lo,hi in [('u_0_60_capped.npz',0,60),('u_26_30_capped.npz',26,30)]:
    p=policy(name)
    assert str(p['model_id'])=='capped_b4_balance'
    assert np.isclose(p['accel_cap'],published_accel_cap()) and np.isclose(p['brake_a'],4)
    assert np.isclose(float(p['t'][-1]),hi)
    assert np.min(p['v'])>=-2e-6 and np.min(p['d'])>=-2e-6
    assert abs(integrate_cells(p,lambda s:feedback_realized(p,s),lo,hi)-float(p['cost']))<2e-4
results=calculate(write=False)
for name in ('U_0_60','U_26_30'):
    r=results[name]
    assert r['value_information_CZK']>0
    assert np.isclose(r['integrated_regret_CZK'],r['feedback_independent_integration_CZK']-r['perfect_information_expected_CZK'],atol=2e-6)
    assert r['regret_range_CZK'][0] > -5e-4
assert not np.isclose(results['U_0_60']['perfect_information_expected_CZK'],4.27575158743,atol=1e-6)
print('Session 5A tests PASS: 66 capped knots, model isolation, residuals, support and regret')
