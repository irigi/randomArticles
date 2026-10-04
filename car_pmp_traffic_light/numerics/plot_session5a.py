"""Checkpoint diagnostic figures; not publication figures."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from session5a_voi import load_known, policy, feedback_realized, MODEL
from value_of_information import load_known as original_load_known
root=Path(__file__).resolve().parents[1]
root.joinpath('diagnostics').mkdir(exist_ok=True)
s,known_cost,kf=load_known()
so,orig_cost,of=original_load_known()
x=np.linspace(s[0],60,480)
fig,ax=plt.subplots(figsize=(8.1,4.6))
ax.plot(x,kf(x),color='#BF6A36',lw=2.1,label='Capped acceleration, braking 4 m/s²')
ax.plot(x,of(x),color='#526E91',lw=1.7,ls='--',label='Original mechanics, braking 7 m/s²')
ax.set(xlabel='Known green time S (s)',ylabel='Best numerical cost (CZK)',title='Known-switch local-best cost comparison')
ax.grid(alpha=.2);ax.legend(frameon=False);fig.tight_layout()
fig.savefig(root/'diagnostics'/'session5a_known_comparison.png',dpi=170)
plt.close(fig)
with np.load(root/'data'/'voi_capped_regret.npz',allow_pickle=False) as z:
    p={k:z[k] for k in z.files}
fig,ax=plt.subplots(figsize=(8.1,4.6))
ax.plot(p['S'],p['regret_broad'],color='#BF6A36',lw=2,label='Broad U(0,60)')
ax.plot(p['late_S'],p['late_regret'],color='#148267',lw=1.6,label='Late U(26,30)')
ax.axhline(0,color='black',lw=.7)
ax.set(xlabel='Realized switch time S (s)',ylabel='Regret (CZK)',title='Capped-model cost of missing perfect switch-time information')
ax.grid(alpha=.2);ax.legend(frameon=False)
fig.tight_layout();fig.savefig(root/'diagnostics'/'session5a_regret.png',dpi=170)
plt.close(fig)
print('Saved capped-known comparison and capped-regret diagnostic PNGs')
