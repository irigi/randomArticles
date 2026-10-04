"""Fast regression and provenance checks; no expensive NLP solve."""
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import quad
HERE=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(HERE/'numerics'))
from model import Model, published_accel_cap, road_load, VMAX
from solve_policies import solve_uniform
D=HERE/'data'

def cost(case, restriction, n):
    a=np.load(D/f's4_{case}_{restriction}_N{n}.npz')
    assert a['t'].size==n+1
    assert a['u'].size==n
    assert np.min(a['v'])>=-2e-5 and np.min(a['d'])>=-2e-5
    if restriction=='before27':
        model=Model(300 if '300' in case else 900,
                    accel_cap=published_accel_cap() if ('capped' in case or 'caponly' in case) else None,
                    brake_a=4 if 'capped' in case else 7)
        t=a['t'][:-1]; selected=t < 27-1e-9
        gap=a['u'][selected]-np.asarray(model.cruise_throttle(a['v'][:-1][selected]))
        assert max(gap) < 2e-5, (case, float(max(gap)))
        assert abs(np.interp(27.,a['t'],a['v'])-a['v'][selected.sum()])<1e-9
        assert max(np.diff(a['v'])[selected]) < 2e-5
    return float(a['cost'])

# The first cell at which acceleration begins is *partial*, not full throttle.
a=np.load(D/'narrow900_N1100.npz')
k=round(26.625/(27.5/1100))
assert abs(a['t'][k]-26.625)<1e-9
assert .65<float(a['u'][k])<.70
assert a['u'][k]<.95 and a['u'][k+1]>.95

for case,n in [('original900',1100),('original300',1100),('capped900',1100),('capped300',1100)]:
    unrestricted=cost(case,'unrestricted',n) if case.startswith('capped') else float(np.load(D/('narrow900_N1100.npz' if case=='original900' else 'narrow300_N1100.npz'))['cost'])
    precut=cost(case,'before27',n)
    assert precut>unrestricted+1e-5,(case,precut,unrestricted)
    if case.endswith('900'):
        all_red=cost(case,'all_red',n)
        assert all_red>precut+1e-5,(case,all_red,precut)
    print(case, 'pre-green advantage',round(precut-unrestricted,9))

for case in ['original900','capped900']:
    hi=cost(case,'before27',1100)
    lo=cost(case,'before27',550)
    assert abs(hi-lo)<4e-5,(case,lo,hi)

baseline=Model(300)
capped=Model(300,accel_cap=published_accel_cap(),brake_a=4)
assert abs(float(baseline.net_full_accel(0))-65.88228)<1e-9
assert abs(float(capped.net_full_accel(0))-published_accel_cap())<1e-9
assert abs(quad(lambda v:1/float(capped.net_full_accel(v)),0,100/3.6)[0]-8.4)<1e-8

# Discretization cutoffs must be aligned and restriction alternatives exclusive.
for kw in [dict(no_accel_before=27.,N=111),dict(no_accel_before=27.,no_speed_increasing_throttle=True,N=1100)]:
    try:
        solve_uniform(27.,27.5,model=baseline,**kw)
    except ValueError:pass
    else: raise AssertionError('solver failed cutoff guard')
print('Session 4: constraints, grid guards, costs, calibration, timing PASS')
