"""Rotation and radial-stability figures, plus surface_results.json.

The other article figures are made by figures.py.
"""
import json,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from drag_regular import make
from drag_second import build
from technical_family import center_data

from paths import DATA, FIG
plt.rcParams.update({'font.family':'serif','font.size':8,'axes.labelsize':9,
 'legend.fontsize':7,'axes.spines.top':False,'axes.spines.right':False,
 'pdf.fonttype':42,'mathtext.fontset':'dejavuserif','lines.linewidth':1.3})
colors=['#125b8a','#b9502b'];styles=['-','--']
def save(fig,name):
 fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',pad_inches=.035)
 fig.savefig(FIG/(name+'.png'),bbox_inches='tight',pad_inches=.035,dpi=180)
 plt.close(fig)

models={a:build(1e-8,a,False,True) for a in [.1,.3]}

threshold=json.loads((DATA/'thresholds.json').read_text())
ar=threshold['monotonicity_candidate']['alpha']

fig,ax=plt.subplots(2,2,figsize=(7.05,4.45),layout='constrained')
for i,a in enumerate(models):
 rot,metric,local,detail=models[a];r=np.linspace(.001,1,700);v=np.array([local(x) for x in r])
 ax[0,0].plot(r,v[:,11],c=colors[i],label=rf'$h_0,\ \alpha={a}$')
 ax[0,0].plot(r,r**3*v[:,6],'--',c=colors[i],label=rf'$m_0/R,\ \alpha={a}$')
 ax[0,1].plot(r,v[:,8],c=colors[i],label=rf'$h_2,\ \alpha={a}$')
 ax[0,1].plot(r,v[:,9]-v[:,8],'--',c=colors[i],label=rf'$k_2,\ \alpha={a}$')
 theta=np.linspace(0,2*np.pi,600);P2=(3*np.cos(theta)**2-1)/2
 epsrot=.05*rot['M'];d0=epsrot*rot['xi0'];d2=epsrot*rot['shape']
 rho=(1+d0+d2*P2)*np.sin(theta)
 zz=(1+d0-d2/2)*np.cos(theta)+1.5*d2*np.cos(theta)**3
 ax[1,0].plot(rho,zz,c=colors[i],label=rf'$\alpha={a}$')
 t=np.linspace(0,np.pi/2,300);p2=(3*np.cos(t)**2-1)/2
 ax[1,1].plot(t*180/np.pi,rot['xi0']+rot['xi2']*p2,c=colors[i],label=rf'$\alpha={a}$')
tt=np.linspace(0,2*np.pi,400);ax[1,0].plot(np.sin(tt),np.cos(tt),':',c='.4',label='Nonrotating')
for p in ax[0]:p.set(xlabel=r'$r/R$',ylabel=r'Coefficient of $(\Omega R)^2$');p.legend(frameon=False)
ax[1,0].set(xlabel=r'$\varpi_{\rm emb}/R$',ylabel=r'$Z_{\rm emb}/R$',aspect='equal');ax[1,0].legend(frameon=False,loc='center',fontsize=6)
ax[1,1].set(xlabel=r'Colatitude $\theta$ (degrees)',ylabel=r'$(r_{\rm surf}/R-1)/(\Omega R)^2$');ax[1,1].legend(frameon=False)
save(fig,'rotation')

seq=list(csv.DictReader(open(DATA/'fixed_eos_sequences.csv')))
rad=json.loads((DATA/'radial_modes.json').read_text());rt=json.loads((DATA/'radial_threshold.json').read_text())
fig,ax=plt.subplots(1,2,figsize=(7.05,2.45),layout='constrained')
for i,a in enumerate([.1,.3,1.,8.]):
 q=[x for x in seq if float(x['alpha'])==a];ref=float(q[-1]['M'])
 ax[0].plot([float(x['hc_fraction']) for x in q],[float(x['M'])/ref for x in q],label=rf'$\alpha_{{\rm ref}}={a:g}$')
ax[0].set(xlabel=r'$h_c/h_{c,\rm ref}$',ylabel=r'$M/M_{\rm ref}$');ax[0].legend(frameon=False)
rr=[x for x in rad if x['degree']==24]+[rt['radial']];rr=sorted(rr,key=lambda x:x['alpha'])
ax[1].semilogx([x['alpha'] for x in rr],[x['sigma2_R2'][0] for x in rr],'o-',ms=3,c=colors[0])
ax[1].axhline(0,ls=':',c='.4');ax[1].axvline(ar,ls='--',c='#177a64',lw=.8)
ax[1].set(xlabel=r'$\alpha$',ylabel=r'$\sigma_0^2 R^2$')
save(fig,'stability')

shape={str(a):{k:float(models[a][0][k]) for k in ['xi0','xi2','shape','k2surface','deltaM','Q','I']} for a in models}
(DATA/'surface_results.json').write_text(json.dumps(shape,indent=2))
