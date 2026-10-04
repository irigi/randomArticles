"""Bounded regression checks for the uniform-prior solver (Session 2).

These checks exercise distinct basins at coarse resolution. Publication-resolution
candidate details live in diagnostics/session2_uncertain_solver.json.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from model import Model
from solve_policies import solve_uniform, _rk4_step

ROOT=Path(__file__).resolve().parents[1]


def check_candidate(result, model):
    t=np.asarray(result['t']); v=np.asarray(result['v']); d=np.asarray(result['d'])
    u=np.asarray(result['u']); b=np.asarray(result['b']); dt=t[-1]/len(u)
    assert np.min(v) >= -1e-6 and np.min(d) >= -1e-6
    defect=0.; minstage=float(np.min(v))
    for k in range(len(u)):
        vn,dn,st=_rk4_step(model,v[k],d[k],u[k],b[k],dt,return_stages=True)
        minstage=min(minstage,*map(float,st))
        defect=max(defect,abs(float(vn)-v[k+1]),abs(float(dn)-d[k+1]))
    assert minstage >= -1e-6 and defect < 2e-5


def main():
    m=Model(300,rest_mode='balance')
    cases=[
        ('broad',0.,60.,240,4.4080),
        ('late',26.,30.,240,3.8432),
        ('narrow',27.,27.5,300,3.6143),
    ]
    report={}
    for name,lo,hi,N,bound in cases:
        r=solve_uniform(lo,hi,N=N,model=m,max_cpu_time=8.0)
        check_candidate(r,m)
        assert r['cost'] < bound
        success=[c for c in r['candidates'] if c.get('success')]
        assert len(success) >= 4
        report[name]={
            'lo':lo,'hi':hi,'N':N,'cost_czk':r['cost'],
            'min_speed_kmh':float(np.min(r['v'])*3.6),
            'selected_start':r['selected_start'],'candidates':r['candidates']}
        print(name,r['cost'],r['selected_start'],len(success),'successful starts',flush=True)
    # The smooth rest convention is not literal pedal behavior. Under the
    # calibrated fuel law gamma(0)=0, so its balancing throttle adds no fuel
    # above the idle term at exactly zero speed.
    report['rest_convention']={
        'balance_unforced_dvdt':float(m.rhs_v(0,0,0)),
        'balance_throttle_at_rest':float(m.cruise_throttle(0.0)),
        'gamma_at_rest':float(m.gamma(0.0)),
    }
    assert abs(report['rest_convention']['gamma_at_rest']) < 1e-14
    out=ROOT/'diagnostics'/'session2_uncertain_coarse_regression.json'
    out.write_text(json.dumps(report,indent=2))
    print('PASS',out)

if __name__=='__main__': main()
