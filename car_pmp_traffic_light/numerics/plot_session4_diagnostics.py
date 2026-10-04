"""Plot verified Session 4 red-branch comparisons (diagnostic, not publication)."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent
D=ROOT/'data'; OUT=ROOT/'diagnostics'
BROWN='#BF6A36'; GRAY='#707984'; CHARCOAL='#303641'
for case in ('original900','capped900','capped300'):
    n=1100
    p={r:np.load(D/f's4_{case}_{r}_N{n}.npz') for r in ['unrestricted','before27'] if (D/f's4_{case}_{r}_N{n}.npz').exists()}
    if case=='original900': p['unrestricted']=np.load(D/'narrow900_N1100.npz')
    if case!='capped300':
        p['all_red']=np.load(D/f's4_{case}_all_red_N{n}.npz')
    fig,ax=plt.subplots(2,1,figsize=(9.2,6.8),sharex=True,gridspec_kw={'height_ratios':[1.35,1]})
    for restriction, label, color, ls, z in (
            ('unrestricted','Unrestricted conditional-red trajectory',BROWN,'-',3),
            ('before27','No acceleration before earliest green',CHARCOAL,'--',2),
            ('all_red','No acceleration throughout red',GRAY,':',1)):
        if restriction not in p: continue
        a=p[restriction];t=a['t']
        ax[0].plot(t,a['v']*3.6,ls=ls,color=color,lw=2.3,zorder=z,label=label)
        ax[1].step(t,np.r_[a['u'],a['u'][-1]],where='post',ls=ls,color=color,lw=1.9,zorder=z)
    for a in ax:
        a.axvspan(27.,27.5,color='#d8e7e2',alpha=.35,zorder=0)
        a.axvline(27.,color=GRAY,lw=1.0,ls=':')
        a.axvline(27.5,color=GRAY,lw=1.0,ls=':')
        a.grid(alpha=.14);a.spines[['top','right']].set_visible(False)
    ax[0].legend(fontsize=8.4,loc='upper left')
    ax[0].set_ylabel('Speed (km/h)');ax[1].set_ylabel('Throttle (0–1)')
    ax[1].set_xlabel('Time since observing red (s)')
    ax[1].set_ylim(-.04,1.04);ax[1].set_xlim(25.8,27.5)
    fig.suptitle(case+' — pre-earliest-green acceleration test (not publication figure)',fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT/f'session4_{case}_comparison.png',dpi=190,bbox_inches='tight')
    plt.close(fig)
print('Wrote 3 Session 4 diagnostic figures')
