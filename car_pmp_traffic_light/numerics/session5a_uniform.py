"""Session 5A: capped broad/late policies, bounded one per CLI invocation."""
import argparse, json
from pathlib import Path
import numpy as np
from model import Model,published_accel_cap
from solve_policies import solve_uniform,save_npz
P=argparse.ArgumentParser()
P.add_argument('--case',choices=['broad','late'],required=True)
P.add_argument('--cpu-per-start',type=float,default=10.)
P.add_argument('--starts',nargs='+',default=['legacy','rolling_0.65','rolling_1.05','coast_then_stop','brake_recover'])
a=P.parse_args()
settings={'broad':(0.,60.,480,'u_0_60_capped.npz'),
          'late':(26.,30.,360,'u_26_30_capped.npz')}
lo,hi,N,out=settings[a.case]
model=Model(300,accel_cap=published_accel_cap(),brake_a=4,rest_mode='balance')
r=solve_uniform(lo,hi,N=N,model=model,start_labels=a.starts,max_cpu_time=a.cpu_per_start)
root=Path(__file__).resolve().parents[1]
save_npz(root/'data'/out,r)
# Preserve the conventional NPZ fields and add an unambiguous model tag.
with np.load(root/'data'/out,allow_pickle=False) as z:
    payload={k:z[k] for k in z.files}
np.savez_compressed(root/'data'/out, **payload, model_id='capped_b4_balance',
    accel_cap=model.accel_cap,brake_a=model.brake_a)
info={'case':a.case,'model':'capped_b4_balance','N':N,'lo':lo,'hi':hi,
      'cost_CZK':float(r['cost']),'min_speed_kmh':float(min(r['v'])*3.6),
      'selected_start':r.get('selected_start'),'candidates':r.get('candidates',[])}
(root/'diagnostics'/f'session5a_{a.case}.json').write_text(json.dumps(info,indent=2)+'\n')
print(json.dumps(info,indent=2),flush=True)
