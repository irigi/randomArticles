"""Fast, no-NLP validation of Session 5B figure inputs and outputs."""
from pathlib import Path
import json
import numpy as np
from PIL import Image
import fitz

ROOT=Path(__file__).resolve().parents[1]
with open(ROOT/'figures'/'session5b_manifest.json') as f:m=json.load(f)
assert m['model_id']=='capped_b4_balance'
assert len(m['files'])==4
with np.load(ROOT/'data'/'u_0_60_capped.npz') as z:
    assert str(z['model_id'])==m['model_id']
    assert abs(float(z['cost'])-m['broad_feedback_expected_cost_CZK'])<1e-10
    assert np.min(z['v'])*3.6<.001
with np.load(ROOT/'data'/'u_26_30_capped.npz') as z:
    assert str(z['model_id'])==m['model_id']
    assert abs(float(z['cost'])-m['late_feedback_expected_cost_CZK'])<1e-10
    assert np.min(z['v'])*3.6>25
with np.load(ROOT/'data'/'s4_capped900_unrestricted_N1100.npz') as z:
    t=z['t'];v=z['v'];u=z['u']; assert len(u)==1100
    assert abs(float(z['cost'])-m['anticipation_cost_unrestricted_CZK'])<1e-10
with np.load(ROOT/'data'/'s4_capped900_before27_N1100.npz') as z:
    assert abs(float(z['cost'])-m['anticipation_cost_no_pre27_acceleration_CZK'])<1e-10
    assert abs(m['anticipation_cost_no_pre27_acceleration_CZK']-m['anticipation_cost_unrestricted_CZK']-m['anticipation_pre27_benefit_CZK'])<1e-11
assert 26.4<=m['anticipation_pre27_acceleration_onset_s']<27.0
with open(ROOT/'diagnostics'/'session5a_voi.json') as f:r=json.load(f)
assert abs(m['broad_VOI_CZK']-r['U_0_60']['value_information_CZK'])<1e-10
assert abs(m['late_VOI_CZK']-r['U_26_30']['value_information_CZK'])<1e-10
for name in m['files']:
    png=ROOT/'figures'/f'{name}.png'
    pdf=ROOT/'figures'/f'{name}.pdf'
    assert png.is_file() and pdf.is_file()
    assert png.stat().st_size>40_000 and pdf.stat().st_size>10_000
    with Image.open(png) as im:
        assert im.width>1200 and im.height>550,(name,im.size)
    with fitz.open(pdf) as doc:
        assert len(doc)==1 and doc[0].rect.width>600
        assert len(doc[0].get_text())>150  # vector/text, not merely a raster
print('PASS: 4 capped publication figure pairs, model provenance, outcomes, onset and VOI consistency')
