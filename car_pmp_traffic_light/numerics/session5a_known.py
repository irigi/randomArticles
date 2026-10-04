"""Session 5A: checkpointed capped-model known-S local-best reference.

Run with e.g. --direction forward or reverse. Every S saves independently;
no old uncapped results are used, mixed, or overwritten.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from model import Model, published_accel_cap, D0, VMAX
from solve_policies import solve_known

P=argparse.ArgumentParser()
P.add_argument('--direction',choices=['forward','reverse'],default='forward')
P.add_argument('--cells',type=int,default=480)
P.add_argument('--cpu-per-start',type=float,default=6.0)
P.add_argument('--times',nargs='*',type=float)
P.add_argument('--extra-starts',action='store_true',help='also challenge with independent legacy/stop starts')
a=P.parse_args()
model=Model(300,accel_cap=published_accel_cap(),brake_a=4,rest_mode='balance')
root=Path(__file__).resolve().parents[1]/'diagnostics'/'capped_known_sweep'
root.mkdir(exist_ok=True,parents=True)
Ts=sorted(set(a.times if a.times else (list(range(22,61))+[26.5,27.5,28.5,29.5])),reverse=(a.direction=='reverse'))
prev=None
for s in Ts:
    fn=root/f'capped_b4_{s:06.2f}_N{a.cells}.npz'
    jp=root/f'capped_b4_{s:06.2f}_N{a.cells}.json'
    old=None
    if fn.exists():
        with np.load(fn,allow_pickle=False) as z:
            old={k:z[k].item() if z[k].ndim==0 else z[k].copy() for k in z.files}
    if old is not None and a.direction=='forward' and not a.extra_starts:
        prev=old
        print(f'S={s:.2f} cached {old["cost"]:.9f}',flush=True)
        continue
    warms=[]
    if prev is not None: warms.append(prev)
    if old is not None: warms.append(old)
    labels=['rolling_0.65','rolling_0.95']
    if a.extra_starts: labels+=['legacy','coast_then_stop','stop_wait']
    labels +=[f'continuation_{j}' for j in range(len(warms))]
    st=time.monotonic()
    try:
        candidate=solve_known(s,N=a.cells,model=model,warm_starts=warms,
                              start_labels=labels,max_cpu_time=a.cpu_per_start)
    except RuntimeError as exc:
        if old is None:
            print(f'S={s:.2f} FAILED: {exc}',flush=True)
            continue
        candidate=None
    if candidate is not None and (old is None or candidate['cost']<old['cost']-1e-9):
        selected=candidate
        np.savez_compressed(fn,**{k:candidate[k] for k in ('t','v','d','u','b','cost','T')},
            model_id='capped_b4_balance',accel_cap=model.accel_cap,brake_a=model.brake_a)
    else: selected=old
    rec={ 'model':'capped_b4_balance','S':s,'cells':a.cells,'direction':a.direction,
          'selected_cost':float(selected['cost']),'selected_min_speed_kmh':float(min(selected['v'])*3.6),
          'trial_cost':float(candidate['cost']) if candidate is not None else None,
          'trial_selected_start':candidate.get('selected_start') if candidate else None,
          'candidates':candidate.get('candidates',[]) if candidate else []}
    if jp.exists():
        earlier=json.loads(jp.read_text());rec['history']=earlier.get('history',[])+[{k:earlier.get(k) for k in ('direction','selected_cost','trial_selected_start','candidates')}]
    jp.write_text(json.dumps(rec,indent=2)+'\n')
    prev=selected
    print(f'S={s:.2f} selected={selected["cost"]:.9f} min_v={min(selected["v"])*3.6:.3f}kmh seconds={time.monotonic()-st:.1f}',flush=True)
