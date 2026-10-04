"""Finite benchmark tests for the deterministic S-known solver (Session 1).

Includes reference comparisons from the independent Astra review. Run:
    python numerics/check_known_multistart.py
The test does NOT certify a global minimizer.
"""
import json
from pathlib import Path
import numpy as np
from model import Model, VMAX, D0
from solve_policies import solve_known, _rk4_step

ROOT=Path(__file__).resolve().parents[1]


def check_candidate(s, result, model):
    v,d,u,b = (np.asarray(result[k]) for k in ('v','d','u','b'))
    dt=s/len(u)
    assert abs(v[0]-VMAX) < 1e-6 and abs(d[0]-D0) < 1e-6
    assert np.min(v) >= -1e-6 and np.max(v) <= VMAX+1e-6
    assert np.min(d) >= -1e-6 and np.max(d) <= D0+1e-6
    assert np.all(np.diff(d) <= 1e-6), 'distance must be nonincreasing'
    max_defect=0
    for k in range(len(u)):
        vnext, dnext, stages = _rk4_step(model, v[k], d[k], u[k], b[k], dt,
                                          return_stages=True)
        assert min(stages) >= -1e-6, 'a negative RK4 intermediate speed'
        max_defect=max(max_defect, abs(vnext-v[k+1]),abs(dnext-d[k+1]))
    assert max_defect < 2e-5
    assert np.isfinite(result['cost'])


def main():
    model=Model(300, rest_mode='balance')
    report={}
    for switch, source_old, bound, floor in [
        (40, 5.060832880168581, 4.9565, 15.),
        (50, 5.947123538336876, 5.9130, 8.)
    ]:
        res=solve_known(switch,N=480,model=model,max_cpu_time=10)
        check_candidate(switch,res,model)
        assert res['cost'] < bound, (switch,res['cost'])
        assert np.min(res['v'])*3.6 > floor
        report[str(switch)]={'old_cost':source_old,
            'new_cost':res['cost'], 'improvement':source_old-res['cost'],
            'min_speed_kmh':float(np.min(res['v'])*3.6),
            'selected_start':res['selected_start'],
            'candidates':res['candidates']}
        print(f"S={switch}: old {source_old:.8f}, new {res['cost']:.8f}; "
              f"min {np.min(res['v'])*3.6:.4f} km/h",flush=True)
    contact=Model(300,rest_mode='contact')
    assert abs(float(contact.rhs_v(0,0,0))) < 1e-12
    assert float(model.rhs_v(0,0,0)) < 0
    assert float(contact.rhs_v(0,1,0)) > 0
    report['rest_mode']={'balance_unforced_dvdt':float(model.rhs_v(0,0,0)),
        'contact_unforced_dvdt':float(contact.rhs_v(0,0,0))}
    out=ROOT/'diagnostics'/'known_solver_regression.json'
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(report,indent=2))
    print(f'PASS: all deterministic solver checks. Report: {out}')

if __name__=='__main__':main()
