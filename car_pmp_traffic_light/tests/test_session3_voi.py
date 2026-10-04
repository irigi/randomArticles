"""Checks complete N=480 reference, support boundaries, and VOI reproducibility.

Run: python tests/test_session3_voi.py
"""
import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'numerics'))
from model import Model,D0,VMAX
from value_of_information import calculate, load_known, load_policy, feedback_realized

knots,cost,fit=load_known()
assert np.all(np.diff(knots)>0)
assert all(np.any(np.isclose(knots,s)) for s in range(22,61))
assert all(np.any(np.isclose(knots,s)) for s in (26.5,27.5,28.5,29.5))
assert abs(fit(21.6)-Model(300).base_route_cost())<1e-10
assert abs(fit(40)-4.9556814696)<1e-6
assert abs(fit(50)-5.9124470932)<1e-6
assert fit(40)<5.06083 and fit(50)<5.94667
result=calculate(write=False)
assert abs(result['U_0_60']['value_information_using_nlp_CZK']-.132081152)<2e-6
assert abs(result['U_26_30']['value_information_using_nlp_CZK']-.171020579)<2e-6
# The repaired calculation must *never* integrate U(26,30) out to S=60.
assert abs(result['U_26_30']['perfect_information_expected_CZK'] - fit.integrate(26,30)/4)<1e-12
assert abs(fit.integrate(26,30)/4 - fit.integrate(26,60)/4)>1
for case in ('U_0_60','U_26_30'):
    assert abs(result[case]['feedback_quadrature_residual_CZK']) < 2e-4
    assert result[case]['minimum_grid_regret_CZK']>-4e-4
with np.load(ROOT/'data'/'voi_regret.npz',allow_pickle=False) as f:
    a={k:f[k] for k in f.files}
assert abs(a['regret_broad'][0])<1e-8
assert 0.4<a['regret_broad'].max()<0.5
assert np.all(a['late_regret']>=0)
assert abs(np.mean(a['S'] == 26)-1/len(a['S']))<1e-6
for fn in ('u_0_60.npz','u_26_30.npz'):
    p=load_policy(fn); assert feedback_realized(p, [p['t'][0],p['t'][-1]]).shape == (2,)
print('PASS: known-time grid, correct prior support, feedback quadrature, regret and VOI')
