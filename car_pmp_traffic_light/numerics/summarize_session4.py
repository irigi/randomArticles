"""Recompute Session 4 metadata from saved NPZ files; no NLP runs."""
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent.parent
D=ROOT/'data'; DIAG=ROOT/'diagnostics'

def load(case,restrict,n=1100):
    if case=='original900' and restrict=='unrestricted': path=D/'narrow900_N1100.npz'
    elif case=='original300' and restrict=='unrestricted': path=D/'narrow300_N1100.npz'
    else: path=D/f's4_{case}_{restrict}_N{n}.npz'
    a=np.load(path)
    return {'cost_czk':float(a['cost']),'file':str(path.relative_to(ROOT))}

summary={}
for case in ('original300','original900','capped300','capped900'):
    row={r:load(case,r) for r in ['unrestricted','before27']}
    if case in ('original300','original900','capped900'):
        row['all_red']=load(case,'all_red')
    u=row['unrestricted']['cost_czk'];pre=row['before27']['cost_czk']
    row['cost_of_disallowing_pre27_acceleration_czk']=pre-u
    if 'all_red' in row:
        row['cost_of_disallowing_all_red_acceleration_czk']=row['all_red']['cost_czk']-u
        row['extra_cost_of_disallowing_acceleration_during_27_to_27p5_czk']=row['all_red']['cost_czk']-pre
    summary[case]=row
# 550-cell cap-only comparison separates cap from brake effect.
caponly={r:load('caponly900',r,550) for r in ('unrestricted','before27')}
caponly['cost_of_disallowing_pre27_acceleration_czk']=caponly['before27']['cost_czk']-caponly['unrestricted']['cost_czk']
summary['caponly900_N550']=caponly
summary['mechanics']={
    'original_net_full_acceleration_at_zero_m_s2':65.88228,
    'capped_net_full_acceleration_at_zero_m_s2':4.136672736266082,
    'original_0_to_100_kmh_s':6.751183406890524,
    'capped_0_to_100_kmh_s':8.4,
    'interpretation':'The capped model is still a simplified acceleration envelope, not a measured powertrain/transient-fuel model.',
}
summary['interpretation']={
  'restriction':'u[k] <= u_cruise(v[k]) only on k with cell start < 27 s; grid exactly includes 27 s. Legacy all_red continues through 27.5 s.',
  'onset_original900_N1100':'First increasing-traction cell starts at 26.625 s with u≈0.676; first near-full-throttle cell starts at 26.650 s.',
  'capped900_N1100':'First increasing-traction and near-full-throttle cell starts at 26.425 s (550-cell onset differs by one discretization cell).',
  'globality':'Best feasible local candidate among tested starts; no global-optimality certificate.',
  'voi':'Session 3 VOI applies ONLY to the original uncapped model. Capped VOI has not been recomputed.',
}
(DIAG/'session4_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
old=json.loads((ROOT/('diagnostics/historical_results_through_session5b.json' if (ROOT/'diagnostics/historical_results_through_session5b.json').exists() else 'results.json')).read_text())
old['narrow_U_27_27_5']['Ct_900'].pop('full_throttle_start_s',None)
old['narrow_U_27_27_5']['Ct_900'].pop('full_throttle_start_position_m',None)
old['narrow_U_27_27_5']['Ct_900']['acceleration_onset_s']=26.625
old['narrow_U_27_27_5']['Ct_900']['throttle_fraction_at_acceleration_onset']=0.6762723
old['narrow_U_27_27_5']['Ct_900']['near_full_throttle_onset_s']=26.65
old['mechanical_sensitivity'].pop('narrow_anticipatory_full_throttle_start_s',None)
old['mechanical_sensitivity']['narrow_anticipatory_near_full_throttle_onset_s_1100_cells']=26.425
old['mechanical_sensitivity']['previous_narrow_anticipation_advantage_is_all_red_restriction']=True
old['session4_anticipation']={k:v for k,v in summary.items() if k!='mechanics' and k!='interpretation'}
old['session4_mechanics']=summary['mechanics']
old['checkpoint_status']['session']=4
old['checkpoint_status']['verified_now'] += ['1100-cell pre-earliest-green and full-red counterfactuals for original model',
  '1100-cell capped model at 300 and 900 CZK/h; 550-cell cap-only comparator',
  'acceleration onset, 0-100 km/h calibration and 550/1100 resolution']
old['checkpoint_status']['pending']=[x for x in old['checkpoint_status']['pending'] if 'pre-earliest' not in x.lower() and 'mechanical' not in x.lower()]
old['checkpoint_status']['pending'] += ['choose final main mechanics, then recompute capped known-S and VOI if promoted',
 'regenerate publication figures and restructure results narrative in Sessions 5-6']
old['cleanup_note']='Session 4 replaces ambiguous acceleration-onset claims, separates before-27s vs all-red costs, and validates capped model; published PDF/figures and overall VOI remain on prior model until Sessions 5-6.'
(ROOT/('diagnostics/historical_results_through_session5b.json' if (ROOT/'diagnostics/historical_results_through_session5b.json').exists() else 'results.json')).write_text(json.dumps(old,indent=2)+'\n')
for case in ('original300','original900','capped300','capped900'):
 print(case,'pre-27 benefit',round(summary[case]['cost_of_disallowing_pre27_acceleration_czk'],9))
