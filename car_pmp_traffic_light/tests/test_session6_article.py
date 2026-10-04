"""Checks current publication article, model tags, expected results and figure paths."""
from pathlib import Path
import json, re
import numpy as np
import fitz
root=Path(__file__).resolve().parents[1]
r=json.loads((root/'results.json').read_text());tex=(root/'article.tex').read_text()
assert r['publication_model']=='capped_b4_balance'
assert (root/'diagnostics/historical_results_through_session5b.json').is_file()
assert '0.184983' in tex and '0.176971' in tex
assert '0.010246' in tex and '0.001385' in tex
assert '4.13309' in tex and '3.84990' in tex
for term in ['65.88','26.425','4.136672736','before the','continuation']:
 assert term in tex,term
for name in ['fig1_known28_capped','fig2_information_capped','fig3_anticipation_capped','fig4_regret_capped']:
 assert f'figures/{name}.pdf' in tex
 assert (root/'figures'/f'{name}.pdf').is_file()
assert not re.search(r'fig\d_[a-z]+_conditional\.png',tex)
assert 'fig2_red_conditional' not in tex
assert abs(r['S28_realized_cost_CZK']['known']-3.6755043875231594)<1e-9
assert abs(r['S28_realized_cost_CZK']['broad_feedback']-4.133089972033455)<1e-8
assert abs(r['S28_realized_cost_CZK']['late_feedback']-3.8499041256993953)<1e-8
for prior,k in [('U_0_60',.18498289999579853),('U_26_30',.17697070417982896)]:
 assert abs(r['capped_value_of_information_CZK'][prior]['value_information_CZK']-k)<1e-10
with fitz.open(root/'article.pdf') as pdf:
 assert len(pdf)==7,len(pdf)
 text=' '.join(page.get_text() for page in pdf)
 assert all(f'Fig. {i}.' in text for i in (1,2,3,4))
 assert all(text.count(f'{s:.6f}')>0 for s in [0.184983,0.176971,0.010246])
 assert pdf[6].get_text().count('References')==1 or 'REFERENCES' in pdf[6].get_text()
print('PASS Session 6: publication PDF, figures, captions, exact numeric outputs, archive semantics')
