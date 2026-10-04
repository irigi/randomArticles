"""Checkpoint-only plots. Paper layout and figure captions belong to Session 5."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parents[1]
D=HERE/'diagnostics'; D.mkdir(exist_ok=True)
with np.load(HERE/'data'/'voi_regret.npz') as z:
    v={k:z[k] for k in z.files}
with np.load(HERE/'data'/'voi_grid.npz') as z:
    x=z['T']; y=z['C']

fig, ax=plt.subplots(figsize=(7.7,3.6))
ax.plot(x,y,'-',color='#404040',lw=1.8,marker='.',markersize=2.5)
ax.axvline(21.6,color='#aaaaaa',ls=':',lw=1.2)
ax.set(xlabel='Realized green-switch time S (s)',ylabel='Cost with known S (CZK)',title='Multistart known-switch reference (480 cells)')
ax.grid(alpha=.18);fig.tight_layout();fig.savefig(D/'session3_known_curve.png',dpi=170);plt.close(fig)

for typ, t, r, title in [
    ('broad',v['S'],v['regret_broad'],'Broad prior U(0, 60)'),
    ('late',v['late_S'],v['late_regret'],'Late prior U(26, 30)'),
]:
    fig,ax=plt.subplots(figsize=(7.7,3.6))
    ax.axhline(0,color='#666666',lw=.8)
    ax.plot(t,r,color='#BF6A36',lw=2.1)
    ax.fill_between(t,0,r,where=r>=0,color='#BF6A36',alpha=.12)
    ax.set(xlabel='Realized green-switch time S (s)',ylabel='Regret R(S) (CZK)',
           title=f'Information regret — {title}')
    ax.grid(alpha=.18);fig.tight_layout()
    fig.savefig(D/f'session3_regret_{typ}.png',dpi=170);plt.close(fig)
