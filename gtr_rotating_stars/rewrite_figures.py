"""Vector figures for the revised two-column manuscript."""
import json,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from drag_regular import make
from drag_second import build
from technical_family import center_data

OUT=Path('output/pdf/revision');FIG=OUT/'figures';FIG.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.size':8,'axes.labelsize':9,
 'legend.fontsize':7,'axes.spines.top':False,'axes.spines.right':False,
 'pdf.fonttype':42,'mathtext.fontset':'dejavuserif','lines.linewidth':1.3})
colors=['#125b8a','#b9502b'];styles=['-','--']
def save(fig,name):
 fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',pad_inches=.035)
 fig.savefig(FIG/(name+'.png'),bbox_inches='tight',pad_inches=.035,dpi=180)
 plt.close(fig)

models={a:build(1e-8,a,False,True) for a in [.1,.3]}
fig,ax=plt.subplots(2,2,figsize=(7.05,4.4),layout='constrained')
for i,a in enumerate(models):
 rot,metric,local,detail=models[a];bg=detail['background'];F,D,A,N,j,aa,fields,fb,sol,c=bg
 r=np.linspace(.0001,.998,1000);s=r*r;f,fp,E,P,v,Es,Ps=fields(s)
 Ec,Pc,*_=center_data(a,c);rho=(E+P)*np.sqrt(f)/(j(np.sqrt(1-s))*np.sqrt(fb));rhoc=(Ec+Pc)/(j(1)*np.sqrt(fb))
 lab=rf'$\alpha={a}$'
 ax[0,0].plot(r,E/Ec,c=colors[i],label=lab+r', $\varepsilon/\varepsilon_c$')
 ax[0,0].plot(r,P/Pc,'--',c=colors[i],label=lab+r', $P/P_c$')
 ax[0,1].plot(rho/rhoc,P/Pc,c=colors[i],label=lab)
 ax[1,0].plot(E/Ec,v,c=colors[i],label=lab)
 ax[1,1].plot(r,(E+P)*v/P,c=colors[i],label=lab)
x=np.linspace(0,1,301)
ax[0,1].plot(x,x**(5/3),':',c='.3',label=r'$x^{5/3}$')
ax[0,1].plot(x,x*x,'-.',c='.55',label=r'$x^2$')
ax[0,0].set(xlabel=r'$r/R$',ylabel='Normalized density and pressure')
ax[0,1].set(xlabel=r'$x=\rho_0/\rho_{0c}$',ylabel=r'$P/P_c$')
ax[1,0].set(xlabel=r'$\varepsilon/\varepsilon_c$',ylabel=r'$c_s^2$')
ax[1,1].set(xlabel=r'$r/R$',ylabel=r'$\Gamma$');ax[1,1].axhline(5/3,ls=':',c='.5')
for p in ax.flat:p.legend(frameon=False)
save(fig,'background')

rows=json.loads(Path('technical_review/family_results.json').read_text())
threshold=json.loads(Path('technical_review/thresholds.json').read_text())
ar=threshold['monotonicity_candidate']['alpha'];ac=threshold['causality_candidate']['alpha'];ast=.3762858476
a=[x['alpha'] for x in rows]
fig,ax=plt.subplots(1,2,figsize=(7.05,2.35),layout='constrained')
ax[0].semilogx(a,[x['C'] for x in rows],'o-',ms=3,c=colors[0]);ax[0].set(xlabel=r'$\alpha$',ylabel=r'$\mathcal{C}=M/R$')
ax[1].semilogx(a,[x['cs2_max'] for x in rows],'o-',ms=3,c=colors[0],label='Maximum')
ax[1].semilogx(a,[x['cs2_center'] for x in rows],'s--',ms=3,c=colors[1],label='Centre')
ax[1].axhline(1,ls=':',c='.4');ax[1].set(xlabel=r'$\alpha$',ylabel=r'$c_s^2$');ax[1].legend(loc='upper left',frameon=False)
for p in ax:
 p.axvspan(.03,ar,alpha=.07,color=colors[0])
 p.axvline(ar,c='#177a64',ls='--',lw=.8);p.axvline(ast,c='#9c5234',ls=':',lw=.8);p.axvline(ac,c='.4',ls='-.',lw=.8)
save(fig,'family')

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

seq=list(csv.DictReader(open('technical_review/fixed_eos_sequences.csv')))
rad=json.loads(Path('technical_review/radial_modes.json').read_text());rt=json.loads(Path('technical_review/radial_threshold.json').read_text())
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

fig,ax=plt.subplots(1,2,figsize=(7.05,2.45),layout='constrained')
use=[x for x in rows if .1<=x['alpha']<=8]
xx=np.logspace(np.log10(1.45),np.log10(20),200)
yy=np.exp(np.polynomial.polynomial.polyval(np.log(xx),[1.35,.697,-.143,.0994,-.0124]))
ax[0].loglog(xx,yy,'--',c='.4',label='Yagi--Yunes fit')
ax[0].loglog([x['Qbar'] for x in use],[x['Ibar'] for x in use],'o',ms=3,c=colors[0],label='Present family')
ax[0].set(xlabel=r'$\bar Q$',ylabel=r'$\bar I$');ax[0].legend(frameon=False)
ax[1].plot([x['C'] for x in use],[100*x['IQ_fractional_residual'] for x in use],'o-',ms=3,c=colors[0],label=r'$I$--$Q$')
ax[1].plot([x['C'] for x in use],[100*x['ILove_fractional_residual'] for x in use],'s--',ms=3,c=colors[1],label=r'$I$--Love')
ax[1].axhline(0,ls=':',c='.5');ax[1].set(xlabel=r'$\mathcal{C}$',ylabel='Fit residual (%)');ax[1].legend(frameon=False)
save(fig,'universal')

shape={str(a):{k:float(models[a][0][k]) for k in ['xi0','xi2','shape','k2surface','deltaM','Q','I']} for a in models}
Path('rewrite/surface_results.json').write_text(json.dumps(shape,indent=2))
