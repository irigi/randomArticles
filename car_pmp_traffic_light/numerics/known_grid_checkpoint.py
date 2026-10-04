"""Resumable known-switch cost sweep, separate from VOI integration.

Example:
    python numerics/known_grid_checkpoint.py --times 25 28 30 35 40 45 50 55 60
    python numerics/known_grid_checkpoint.py --start 22 --end 60 --step 1

Individual switch solutions and diagnostics are saved *after each solve*. Re-run
with --reverse to warm-start from larger times and challenge the current best
with independent multistart candidates. This script is intended for a batched,
checkpointed numerical campaign, not proof of global optimality.
"""
import argparse, json
from pathlib import Path
import numpy as np
from model import Model, VMAX, D0
from solve_policies import solve_known

ap=argparse.ArgumentParser()
ap.add_argument('--times', type=float, nargs='*', default=[])
ap.add_argument('--start', type=float)
ap.add_argument('--end', type=float)
ap.add_argument('--step', type=float, default=1)
ap.add_argument('--reverse', action='store_true')
ap.add_argument('--starts', nargs='+', default=['legacy','rolling_0.65'])
ap.add_argument('--cells', type=int, default=480)
ap.add_argument('--cpu-per-start', type=float, default=10)
ap.add_argument('--rest-mode', choices=['balance','contact'], default='balance')
args=ap.parse_args()
if args.start is not None:
    if args.end is None: ap.error('--end required with --start')
    times=list(np.arange(args.start,args.end+args.step*.25,args.step))
else: times=args.times
if not times: ap.error('specify --times or --start/--end')
times=sorted(set(times), reverse=args.reverse)
model=Model(300,rest_mode=args.rest_mode)
root=Path(__file__).resolve().parents[1]/'diagnostics'/'known_sweep'
root.mkdir(parents=True,exist_ok=True)
prev=None
for tau in times:
    fn=root/f'{args.rest_mode}_{tau:06.2f}_N{args.cells}.npz'
    diag=root/f'{args.rest_mode}_{tau:06.2f}_N{args.cells}.json'
    old=None
    if fn.exists():
        with np.load(fn) as z:
            old={k:z[k].item() if z[k].ndim == 0 else z[k].copy() for k in z.files}
    # The opposite-direction pass challenges stored results; same-direction
    # default skips a completed time, so interrupted batches are resumable.
    if old is not None and not args.reverse:
        prev=old
        print(f'S={tau:.2f}: cached {old["cost"]:.9f}',flush=True)
        continue
    warms=[]
    if prev is not None: warms.append(prev)
    if old is not None: warms.append(old)
    result=solve_known(tau,N=args.cells,model=model,max_cpu_time=args.cpu_per_start,
                       warm_starts=warms,
                       start_labels=args.starts + [f"continuation_{j}" for j in range(len(warms))])
    if old is None or float(result['cost']) < float(old['cost'])-1e-9:
        np.savez(fn,**{k:result[k] for k in ('t','v','d','u','b','cost','T')})
        selected=result
    else: selected=old
    diag.write_text(json.dumps({'S':tau,'cost':float(result['cost']),
        'selected_cost':float(selected['cost']),
        'selected_start':result.get('selected_start'),
        'candidates':result.get('candidates',[])},indent=2))
    prev=selected
    print(f'S={tau:.2f}: selected {float(selected["cost"]):.9f}, '
          f'trial {float(result["cost"]):.9f}; '
          f'min speed {min(selected["v"])*3.6:.4f} km/h',flush=True)
print('Completed checkpointed sweep in',root)
