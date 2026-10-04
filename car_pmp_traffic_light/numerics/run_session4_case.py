"""Independently resumable Session 4: one narrow-prior solve per invocation.

For pre-green tests the earliest-green boundary is exactly a grid node.
Never batch IPOPT starts/cases: each invocation is a checkpoint.
"""
import argparse, json
from pathlib import Path
import numpy as np
from model import Model, published_accel_cap, VMAX
from solve_policies import solve_uniform, save_npz

CASES = {
 'original300': (300,None,7.0),
 'original900': (900,None,7.0),
 'capped300': (300,published_accel_cap(),4.0),
 'capped900': (900,published_accel_cap(),4.0),
 'caponly900': (900,published_accel_cap(),7.0),
}


def onset(policy, model):
    t=np.asarray(policy['t']); u=np.asarray(policy['u']);v=np.asarray(policy['v']);b=np.asarray(policy['b'])
    throttle_eq=np.asarray(model.cruise_throttle(v[:-1]))
    # Cellwise throttle-above-equilibrium onset, not first full throttle.
    active=np.where((u>throttle_eq+1e-4)&(b<1e-3))[0]
    full=np.where(u>=.95)[0]
    return {'traction_above_cruise_first_cell_s':float(t[active[0]]) if len(active) else None,
            'throttle_first_ge_0p95_cell_s':float(t[full[0]]) if len(full) else None,
            'earliest_green_speed_kmh':float(np.interp(27.,t,v)*3.6),
            'earliest_green_distance_m':float(np.interp(27.,t,np.asarray(policy['d']))),
            'min_speed_kmh':float(min(v)*3.6),
            'max_pre_green_speed_change_cell_m_s':float(np.max(np.diff(v)[t[:-1]<27-1e-9]))}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--case',choices=CASES,required=True)
    p.add_argument('--restriction',choices=['unrestricted','before27','all_red'],required=True)
    p.add_argument('--cells',type=int,choices=[275,550,825,1100],default=550)
    p.add_argument('--starts',nargs='+',default=['legacy','brake_recover'])
    p.add_argument('--limit',type=float,default=12.0)
    args=p.parse_args()
    time_cost,cap,brake=CASES[args.case]
    model=Model(time_cost,accel_cap=cap,brake_a=brake)
    if args.restriction=='before27' and abs(27*args.cells/27.5-round(27*args.cells/27.5))>1e-8:
        p.error('chosen cell count does not put 27 s on a grid node')
    r=solve_uniform(27.,27.5,N=args.cells,model=model,
       no_speed_increasing_throttle=args.restriction=='all_red',
       no_accel_before=27. if args.restriction=='before27' else None,
       max_cpu_time=args.limit,start_labels=args.starts)
    out=Path(__file__).resolve().parent.parent
    name=f's4_{args.case}_{args.restriction}_N{args.cells}'
    path=out/'data'/(name+'.npz')
    save_npz(path,r)
    info={'case':args.case,'restriction':args.restriction,'N':args.cells,
          'time_cost_czk_h':time_cost,'accel_cap_m_s2':cap,'brake_a_m_s2':brake,
          'expected_cost_czk':float(r['cost']),'starts':r['candidates'],
          'selected_start':r.get('selected_start'), 'onset':onset(r,model)}
    (out/'diagnostics'/(name+'.json')).write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps({'file':str(path.name),'cost':info['expected_cost_czk'],
          'selected_start':info['selected_start'], 'onset':info['onset']},indent=2),flush=True)

if __name__=='__main__':
    main()
