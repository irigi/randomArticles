"""Fast packaging metadata for verified capped-model Session 5A results.
No NLP and no changes to original-model results or publication figures.
"""
import json
from pathlib import Path
import numpy as np
from session5a_voi import calculate,load_known,ROOT

r=calculate(write=True)
original=json.loads((ROOT/'diagnostics'/'session3_voi.json').read_text())
s,y,known=load_known()
broad,late=r['U_0_60'],r['U_26_30']
assert len(s)==66
sweep=ROOT/'diagnostics'/'capped_known_sweep'
files=list(sweep.glob('capped_b4_*_N480.npz'))
assert len(files)==43, len(files)
reverse_log=(ROOT/'diagnostics'/'session5a_known_reverse.log').read_text()
reverse_count=sum(l.startswith('S=') for l in reverse_log.splitlines())
assert reverse_count>=18,(reverse_count, 'reverse challenge unfinished')

resolution={}
for tau in (60.,):
    sample={}
    for n in (480,960):
        with np.load(sweep/f'capped_b4_{tau:06.2f}_N{n}.npz',allow_pickle=False) as z:
            sample[str(n)]=float(z['cost'])
    sample['difference_960_minus_480_CZK']=sample['960']-sample['480']
    resolution[str(tau)]=sample

orig_results=json.loads((ROOT/('diagnostics/historical_results_through_session5b.json' if (ROOT/'diagnostics/historical_results_through_session5b.json').exists() else 'results.json')).read_text())
orig_results['capped_main_candidate_Ct_300']={
    'status':'verified Session 5A numerical estimate, not published until figures and manuscript are rebuilt',
    'model_id':'capped_b4_balance','net_accel_cap_m_s2':4.136672736266082,
    'brake_parameter_m_s2':4.0,'known_reference_cells':480,
    'feedback_cells_broad':480,'feedback_cells_late':360,
    'known_S_28_CZK':float(known(28)),
    'known_S_40_CZK':float(known(40)),
    'known_S_60_CZK':float(known(60)),
    'U_0_60':broad,'U_26_30':late,
    'S_60_resolution':resolution['60.0'],
}
orig_results['checkpoint_status']['session']='5A'
orig_results['checkpoint_status']['verified_now'] += [
    'Capped-model broad/late multistart policy regeneration with explicit model tagging',
    'Capped-model 480-cell known-S grid, 43 points plus free-flow constant segment',
    'Capped-model reverse challenge at 18 known switch times',
    'Capped-model VOI and realized-cost quadrature checks',
    'Capped known S=60 480 vs 960 cell resolution check']
orig_results['checkpoint_status']['pending']=[p for p in orig_results['checkpoint_status']['pending']
    if p not in ('choose final main mechanics, then recompute capped known-S and VOI if promoted',)]
orig_results['checkpoint_status']['pending']+=['Publication Fig. 1–4 still use old main model; regenerate in Session 5B',
    'Rebuild manuscript text, citations, tables and PDF in Session 6',
    'IPOPT multistart yields best local candidates, not globally certified optima']
(ROOT/('diagnostics/historical_results_through_session5b.json' if (ROOT/'diagnostics/historical_results_through_session5b.json').exists() else 'results.json')).write_text(json.dumps(orig_results,indent=2)+'\n')
summary={'model_id':'capped_b4_balance',
         'comparison_original_uncapped':{
             name:{'perfect_information_expected_CZK':original[name]['perfect_information_expected_CZK'],
                   'feedback_nlp_expected_CZK':original[name]['feedback_nlp_expected_CZK'],
                   'VOI_CZK':original[name]['value_information_using_nlp_CZK']}
             for name in ('U_0_60','U_26_30')},
         'capped':r,'known_S_28_CZK':float(known(28)),
         'known_S_40_CZK':float(known(40)),'known_S_60_CZK':float(known(60)),
         'known_grid_480_points':len(files),'reverse_challenge_points':reverse_count,
         'resolution':resolution}
(ROOT/'diagnostics'/'session5a_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
readme=f'''# Session 5A — capped mechanics, known-time reference and value of information

This is a **numerical checkpoint only**. The publication `article.tex`, `article.pdf`
and four figures remain from the previous, **uncapped** main example. Do not cite their
VOI table or captions as capped-model results. Session 5B will regenerate plots;
Session 6 will revise the manuscript and rebuild the PDF.

## Selected model

The same calibrated **net full-throttle acceleration cap** of 4.136672736 m/s²
from Session 4, 4 m/s² maximum brake parameter, 300 CZK/h time value,
and smooth `rest_mode=balance` convention (not literal pedal behavior).
`Model.H(v)` is recomputed for this dynamics model. No uncapped green cost or
known-time references are reused.

## Capped model (CZK, best feasible *local* candidates)

| Prior | Feedback expected | Perfect info expected | VOI | Independently integrated regret |
|---|---:|---:|---:|---:|
| Broad U(0,60) | {broad['feedback_nlp_expected_CZK']:.9f} | {broad['perfect_information_expected_CZK']:.9f} | **{broad['value_information_CZK']:.9f}** | {broad['integrated_regret_CZK']:.9f} |
| Late U(26,30) | {late['feedback_nlp_expected_CZK']:.9f} | {late['perfect_information_expected_CZK']:.9f} | **{late['value_information_CZK']:.9f}** | {late['integrated_regret_CZK']:.9f} |

Both priors' feedback objective is independently integrated by cellwise Gauss
quadrature; differences from NLP are **{broad['quadrature_residual_CZK']:.3g}** CZK
and **{late['quadrature_residual_CZK']:.3g}** CZK, respectively.

## Original model vs capped mechanics

| Prior | Original VOI | Capped VOI |
|---|---:|---:|
| U(0,60) | {original['U_0_60']['value_information_using_nlp_CZK']:.6f} | {broad['value_information_CZK']:.6f} |
| U(26,30) | {original['U_26_30']['value_information_using_nlp_CZK']:.6f} | {late['value_information_CZK']:.6f} |

The ordering changes: in the capped model, **broad > late**, while with the
original uncapped mechanics, **late > broad**. These two priors also differ
in mean and support, so this is **not** evidence for a universal trend with
uncertainty width. The broad-policy cost is more sensitive to recovery from
stopping when low-speed acceleration is capped.

## Reference, robustness and provenance

- 43 independently checkpointed capped known-`S` candidates, integers 22..60
  and half-second values 26.5,27.5,28.5,29.5 at **N=480**; plus the constant
  free-flow optimum on `[0,21.6]` (66 interpolation knots in total).
- Forward continuation plus an independently challenged **{reverse_count}-point** reverse
  subset, using physically interpretable rolling starts and cached solutions.
- Capped known costs at S=28,40,60 s: **{float(known(28)):.9f}**,
  **{float(known(40)):.9f}**, **{float(known(60)):.9f}** CZK.
- At S=60 s, N=480 yields {resolution['60.0']['480']:.9f} CZK and N=960 yields
  {resolution['60.0']['960']:.9f} CZK, a difference of
  **{resolution['60.0']['difference_960_minus_480_CZK']:.7f} CZK**.
  This is a spot check, not a dense convergence or globality certificate.
- `data/u_0_60_capped.npz`, `data/u_26_30_capped.npz` and their solver diagnostics
  are explicitly model tagged. Existing uncapped files, original-session VOI
  plots and tables are preserved unchanged for historical reproduction.
- `data/voi_capped_grid.npz`, `data/voi_capped_regret.npz` and
  `diagnostics/session5a_voi.json` are model-specific. No stale references used.

### Reproduction

```bash
OPENBLAS_NUM_THREADS=1 python numerics/session5a_uniform.py --case broad
OPENBLAS_NUM_THREADS=1 python numerics/session5a_uniform.py --case late
OPENBLAS_NUM_THREADS=1 python numerics/session5a_known.py --direction forward
OPENBLAS_NUM_THREADS=1 python numerics/session5a_known.py --direction reverse --times 22 25 26 26.5 27 27.5 28 28.5 29 29.5 30 33 35 40 45 50 55 60
OPENBLAS_NUM_THREADS=1 python numerics/session5a_known.py --times 60 --cells 960
python numerics/session5a_voi.py
python numerics/plot_session5a.py
python tests/test_session5a_capped_voi.py
```

**Limitations:** these costs are best found IPOPT *local* candidates, not
certified optima. The 480-cell reference, interpolation between switch times,
and stage dynamics introduce numerical error, though independent quadrature,
multistart, reverse challenges and the 960-cell endpoint test bound selected
sources. No fixed-mean uncertainty-width sweep has been performed.
'''
(ROOT/'SESSION5A.md').write_text(readme)
readme_file=ROOT/'README.md'
text=readme_file.read_text()
if '## Session 5A checkpoint' not in text:
    readme_file.write_text('## Session 5A checkpoint\n\n**Capped-model numerical results** are in `SESSION5A.md`, '
      '`diagnostics/session5a_summary.json` and `data/voi_capped_*.npz`. '
      'The PDF and numbered publication figures have deliberately **not** been updated yet.\n\n'+text)
print(json.dumps(summary,indent=2))
