"""Stage 2 vector figures, designed for 181-mm two-column width."""
import json,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'review2/stage2'
DATA=BASE/'data';FIG=ROOT/'output/pdf/final/figures';FIG.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.size':8,'axes.labelsize':9,
 'legend.fontsize':7,'axes.spines.top':False,'axes.spines.right':False,
 'pdf.fonttype':42,'mathtext.fontset':'dejavuserif','lines.linewidth':1.4})
COLORS={.1:'#125b8a',.3:'#b9502b'}

def save(fig,name):
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',pad_inches=.065)
    fig.savefig(FIG/(name+'.png'),bbox_inches='tight',pad_inches=.065,dpi=180)
    plt.close(fig)

def load_rows(name):
    with (DATA/name).open() as f:
        return [{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]

def profiles():
    rows=load_rows('thermodynamic_profiles.csv')
    summaries=json.loads((DATA/'thermodynamic_checks.json').read_text())
    sets={a:[v for v in rows if v['alpha']==a] for a in COLORS}
    centers={a:next(v['center'] for v in summaries if v['alpha']==a) for a in COLORS}
    fig,ax=plt.subplots(2,2,figsize=(7.05,4.5),layout='constrained')
    for a,v in sets.items():
        r=np.array([q['r'] for q in v]);center=centers[a];c=COLORS[a];lab=rf'$\alpha={a}$'
        ax[0,0].plot(r,[q['E']/center['E'] for q in v],c=c,label=lab+r', $\varepsilon/\varepsilon_c$')
        ax[0,0].plot(r,[q['P']/center['P'] for q in v],'--',c=c,label=lab+r', $P/P_c$')
        ax[0,1].plot([q['rho0']/center['rho0'] for q in v],[q['P']/center['P'] for q in v],c=c,label=lab)
        ax[1,0].plot([q['E']/center['E'] for q in v],[q['cs2'] for q in v],c=c,label=lab)
        ax[1,1].plot(r,[q['Gamma'] for q in v],c=c,label=lab)
    x=np.linspace(0,1,301)
    ax[0,1].plot(x,x**(5/3),':',c='.3',label=r'$x^{5/3}$')
    ax[0,1].plot(x,x*x,'-.',c='.55',label=r'$x^2$')
    ax[0,0].set(xlabel=r'$r/R$',ylabel=r'$\varepsilon/\varepsilon_c,\ P/P_c$')
    ax[0,1].set(xlabel=r'$x=\rho_0/\rho_{0c}$',ylabel=r'$P/P_c$')
    ax[1,0].set(xlabel=r'$\varepsilon/\varepsilon_c$',ylabel=r'$c_s^2$')
    ax[1,1].set(xlabel=r'$r/R$',ylabel=r'$\Gamma$')
    ax[1,1].axhline(5/3,ls=':',c='.5',label=r'$5/3$')
    ax[1,1].plot([1],[5/3],'o',mfc='white',mec='.2',ms=3,zorder=5)
    for p in ax.flat:p.legend(frameon=False);p.grid(alpha=.12)
    # A small surface inset makes the true approach visible without flattening it.
    ins=ax[1,1].inset_axes([.09,.43,.40,.34])
    for a,v in sets.items():
        vv=[q for q in v if q['r']>=.98]
        ins.plot([q['r'] for q in vv],[q['Gamma'] for q in vv],c=COLORS[a],lw=1)
    ins.axhline(5/3,ls=':',c='.5',lw=.7);ins.set_xlim(.98,1.0002);ins.set_ylim(1.638,1.674)
    ins.set_xticks([.98,1]);ins.set_yticks([1.64,5/3],['1.64',r'$5/3$']);ins.tick_params(labelsize=6)
    ins.spines['top'].set_visible(True);ins.spines['right'].set_visible(True)
    ins.set_title('Surface',fontsize=6,pad=2)
    save(fig,'background')

    fig,ax=plt.subplots(2,2,figsize=(7.05,4.55),layout='constrained')
    for a,v in sets.items():
        r=[q['r'] for q in v];lab=rf'$\alpha={a}$';c=COLORS[a]
        ax[0,0].plot(r,[q['rho0']/centers[a]['rho0'] for q in v],c=c,label=lab)
        ax[0,1].plot(r,[q['enthalpy'] for q in v],c=c,label=lab)
        ax[1,0].plot(r,[q['internal'] for q in v],c=c,label=lab)
        ax[1,1].plot(r,[q['K'] for q in v],c=c,label=lab)
    ax[0,0].set(ylabel=r'$\rho_0/\rho_{0c}$',title='Rest-mass density')
    ax[0,1].set(ylabel=r'$\mathcal{H}_b/(m_b c^2)$',title='Enthalpy per particle')
    ax[1,0].set(ylabel=r'$(\varepsilon-\rho_0 c^2)/(\rho_0 c^2)$',title='Internal energy above rest mass')
    ax[1,1].set(ylabel=r'$K_{\rm eff}/K_{{\rm eff},c}$',title=r'Effective coefficient: $P/\rho_0^{5/3}$')
    for p in ax.flat:p.set(xlabel=r'$r/R$',xlim=(0,1));p.legend(frameon=False);p.grid(alpha=.12)
    ax[0,0].set_ylim(bottom=0);ax[1,0].set_ylim(bottom=0)
    save(fig,'thermodynamic_profiles')

    fig,ax=plt.subplots(1,2,figsize=(7.05,2.4),layout='constrained')
    for a,v in sets.items():
        r=[q['r'] for q in v];lab=rf'$\alpha={a}$';c=COLORS[a]
        ax[0].plot(r,[q['theta']/centers[a]['theta'] for q in v],c=c,label=lab)
        ax[1].plot(r,[q['theta'] for q in v],c=c,label=lab)
    ax[0].set(ylabel=r'$\Theta/\Theta_c$',title='Normalized pressure per particle')
    ax[1].set(ylabel=r'$\Theta=P/(m_b n c^2)$',title='Pressure per particle in rest-energy units')
    for p in ax:p.set(xlabel=r'$r/R$',xlim=(0,1),ylim=(0,None));p.legend(frameon=False);p.grid(alpha=.12)
    save(fig,'pressure_per_particle')

    fig,ax=plt.subplots(1,2,figsize=(7.05,2.65),layout='constrained')
    for a,v in sets.items():
        outer=[q for q in v if q['r']>.975];near=[q for q in v if q['z']<.035]
        ax[0].plot([q['r'] for q in outer],[q['Gamma'] for q in outer],c=COLORS[a],label=rf'$\alpha={a}$')
        ax[1].plot([q['z'] for q in near],[q['Gamma'] for q in near],c=COLORS[a],label=rf'$\alpha={a}$')
    z=np.linspace(0,.035,100);ax[1].plot(z,5/3-5*z/9,':',c='.2',label=r'$5/3-5z/9$')
    ax[0].axhline(5/3,ls=':',c='.4',label=r'$5/3$')
    ax[0].set(xlabel=r'$r/R$',ylabel=r'$\Gamma$',xlim=(.975,1.0003),title='Outer local minimum')
    ax[1].set(xlabel=r'$z=\sqrt{1-r^2/R^2}$',ylabel=r'$\Gamma$',xlim=(0,.035),title='Surface asymptote')
    for p in ax:p.legend(frameon=False);p.grid(alpha=.12)
    save(fig,'gamma_surface_validation')

def universality():
    rows=json.loads((ROOT/'technical_review/family_results.json').read_text())
    refs=json.loads((ROOT/'review2/reference_model_residuals.json').read_text())
    for v in refs:
        v['Lambda']=2*v['k2']/(3*v['C']**5)
        fit=np.exp(np.polynomial.polynomial.polyval(np.log(v['Lambda']),[1.47,.0817,.0149,.000287,-.0000364]))
        v['ILove_fractional_residual']=v['Ibar']/fit-1
    (DATA/'reference_model_residuals.json').write_text(json.dumps(refs,indent=2))
    use=[v for v in rows if .1<=v['alpha']<=8]
    fig,ax=plt.subplots(1,3,figsize=(7.05,2.65),layout='constrained')
    xx=np.geomspace(1.4,22,200);yy=np.exp(np.polynomial.polynomial.polyval(np.log(xx),[1.35,.697,-.143,.0994,-.0124]))
    ax[0].loglog(xx,yy,'--',c='.5',lw=1,label='Yagi–Yunes fit')
    mappings=[('Qbar','Ibar'),('C','IQ_fractional_residual'),('C','ILove_fractional_residual')]
    for i,(xkey,ykey) in enumerate(mappings):
        factor=1 if i==0 else 100
        ax[i].plot([v[xkey] for v in use],[factor*v[ykey] for v in use],'-',c=COLORS[.1],lw=1,zorder=2)
        for v in use:
            ax[i].plot(v[xkey],factor*v[ykey],'o',mec=COLORS[.1],mfc=COLORS[.1] if v['alpha']<.37628585 else 'white',ms=3.5,zorder=4)
        for model,c,mark in [('tolman_vii','#39785c','^'),('polytrope',COLORS[.3],'s'),('buchdahl','#805599','D')]:
            vs=[v for v in refs if v['model']==model]
            ax[i].plot([v[xkey] for v in vs],[factor*v[ykey] for v in vs],linestyle='none',marker=mark,
                       c=c,mfc='white',mew=1,ms=4,zorder=5)
        ax[i].grid(alpha=.12)
    ax[0].set(xlabel=r'$\bar Q$',ylabel=r'$\bar I$');ax[0].set_yticks([5,10,20,40],['5','10','20','40'])
    ax[0].set_xticks([2,5,10,20],['2','5','10','20'])
    ax[1].set(xlabel=r'$\mathcal{C}=M/R$',ylabel=r'$I$–$Q$ residual (%)',ylim=(-6,2.2))
    ax[2].set(xlabel=r'$\mathcal{C}=M/R$',ylabel=r'$I$–Love residual (%)',ylim=(-3.5,2.2))
    for p in ax[1:]:p.axhline(0,ls=':',c='.5',lw=.8);p.set_xticks([.1,.2,.27])
    handles=[Line2D([],[],color=COLORS[.1],marker='o',ms=3,label='Present family'),
             Line2D([],[],color='#39785c',marker='^',mfc='white',linestyle='',label='Tolman VII'),
             Line2D([],[],color=COLORS[.3],marker='s',mfc='white',linestyle='',label=r'$n=1$ (rest-mass EOS)'),
             Line2D([],[],color='#805599',marker='D',mfc='white',linestyle='',label='Buchdahl'),
             Line2D([],[],color='.5',linestyle='--',label='Yagi–Yunes fit')]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False,fontsize=7)
    save(fig,'universal')

def benchmarks():
    cm=json.loads((DATA/'benchmark_uniform.json').read_text());yy=json.loads((DATA/'benchmark_yagi_yunes.json').read_text())
    fig,ax=plt.subplots(1,2,figsize=(7.05,2.65),layout='constrained')
    for key,c,marker in [('I',COLORS[.1],'o'),('Q',COLORS[.3],'s')]:
        ax[0].plot([v['C'] for v in cm],[100*v[key+'_relative_error'] for v in cm],marker,c=c,ms=4,label=rf'$\bar {key}$')
        ax[1].errorbar([v['C'] for v in yy],[100*v[key+'_relative_error'] for v in yy],
            yerr=[100*v[key+'_digitization_relative_allowance'] for v in yy],fmt=marker,c=c,ms=4,capsize=3,label=rf'$\bar {key}$')
    ax[0].set(title='Uniform-density table',ylabel='Calculated / published − 1 (%)')
    ax[1].set(title=r'Published $n=1$ curves',ylabel='Calculated / digitized − 1 (%)')
    for p in ax:p.axhline(0,c='.5',ls=':',lw=.8);p.set_xlabel(r'$\mathcal{C}=M/R$');p.legend(frameon=False);p.grid(alpha=.12)
    save(fig,'external_benchmarks')

def family():
    rows=json.loads((ROOT/'technical_review/family_results.json').read_text())
    t=json.loads((ROOT/'technical_review/thresholds.json').read_text());ar=t['monotonicity_candidate']['alpha'];ac=t['causality_candidate']['alpha'];ast=.3762858476
    a=[v['alpha'] for v in rows];fig,ax=plt.subplots(1,2,figsize=(7.05,2.5),layout='constrained')
    ax[0].semilogx(a,[v['C'] for v in rows],'o-',ms=3,c=COLORS[.1]);ax[0].set(ylabel=r'$\mathcal{C}=M/R$')
    ax[1].semilogx(a,[v['cs2_max'] for v in rows],'o-',ms=3,c=COLORS[.1],label='Maximum')
    ax[1].semilogx(a,[v['cs2_center'] for v in rows],'s--',ms=3,c=COLORS[.3],label='Centre')
    ax[1].axhline(1,ls=':',c='.4');ax[1].set(ylabel=r'$c_s^2$');ax[1].legend(frameon=False,loc='upper left')
    for p in ax:
        p.axvspan(.03,ar,alpha=.06,color=COLORS[.1]);p.axvline(ar,c='#177a64',ls='--',lw=.8)
        p.axvline(ast,c='#9c5234',ls=':',lw=.8);p.axvline(ac,c='.4',ls='-.',lw=.8);p.set_xlabel(r'$\alpha$')
    save(fig,'family')

if __name__=='__main__':
    profiles();universality();benchmarks();family()
