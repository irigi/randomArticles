"""Session-2 checkpoint: regenerate robust unrestricted uniform-prior policies.

Each case is independently runnable/savable so a difficult IPOPT start cannot
invalidate the whole checkpoint. This deliberately excludes the red-acceleration
counterfactual and mechanical sensitivity runs, scheduled for Session 4.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from model import Model
from solve_policies import solve_uniform, save_npz

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
OUT=ROOT/'data'; OUT.mkdir(exist_ok=True)
DIAG=ROOT/'diagnostics'; DIAG.mkdir(exist_ok=True)
DIAGFILE=DIAG/'session2_uncertain_solver.json'

CASES = {
    'broad': ('u_0_60.npz', 0.0, 60.0, 480, 300.0,
              ['legacy','rolling_1.05','stop_wait','coast_then_stop','brake_recover']),
    'late': ('u_26_30.npz', 26.0, 30.0, 360, 300.0,
             ['legacy','rolling_0.65','rolling_1.05','stop_wait','coast_then_stop','brake_recover']),
    # Coarse six-start tests show a common basin. At N=1100 these two
    # qualitatively distinct starts are a bounded publication-resolution pair.
    'narrow300': ('narrow300_N1100.npz', 27.0, 27.5, 1100, 300.0,
                  ['legacy','brake_recover']),
    'narrow900': ('narrow900_N1100.npz', 27.0, 27.5, 1100, 900.0,
                  ['legacy','brake_recover']),
}

ap=argparse.ArgumentParser()
ap.add_argument('--case', choices=[*CASES,'all'], default='all')
args=ap.parse_args()
keys=list(CASES) if args.case=='all' else [args.case]
summary=json.loads(DIAGFILE.read_text()) if DIAGFILE.exists() else {}

for key in keys:
    name,lo,hi,N,ct,labels=CASES[key]
    print(f'solving {key}: {name} ...', flush=True)
    r=solve_uniform(lo,hi,N=N,model=Model(ct),start_labels=labels,
                    max_cpu_time=12.0)
    save_npz(OUT/name,r)
    summary[name]={
        'lo':lo,'hi':hi,'N':N,'time_cost_czk_h':ct,
        'cost_czk':float(r['cost']),
        'selected_start':r.get('selected_start'),
        'min_speed_kmh':float(np.min(r['v'])*3.6),
        'end_distance_m':float(r['d'][-1]),
        'end_speed_kmh':float(r['v'][-1]*3.6),
        'candidates':r.get('candidates',[]),
    }
    DIAGFILE.write_text(json.dumps(summary,indent=2))
    print(summary[name], flush=True)
    print('saved',OUT/name,'and',DIAGFILE,flush=True)
