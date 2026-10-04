from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'numerics'))
from model import Model, VMAX
from solve_policies import branch, integrate_red

DATA=ROOT/'data'; FIG=ROOT/'figures'; FIG.mkdir(exist_ok=True)
GRAY='#303641'; BROWN='#BF6A36'; BROAD=BROWN; LATE=BROWN; NARROW=BROWN; COUNTER='#6E737B'; LIGHT='#A23B3B'; GREEN='#148267'

def load(name,T=None,cost_key='cost'):
    a=np.load(DATA/name)
    r={k:a[k] for k in a.files}
    if T is not None:r['T']=T
    elif 'T' not in r:r['T']=float(r['t'][-1])
    if cost_key in r:r['cost']=float(r[cost_key])
    elif 'cost' in r:r['cost']=float(r['cost'])
    return r

base=Model(300)
known=load('k_28.npz',28)
broad=load('u_0_60.npz',60)
late=load('u_26_30.npz',30)
kb=branch(known,28,base); bb=branch(broad,28,base); lb=branch(late,28,base)

# Fig. 1: two priors, same realized switch, with time panels so waiting is visible.
fig,axs=plt.subplots(2,2,figsize=(12.2,7.2),sharex='row',gridspec_kw={'height_ratios':[1.25,1.0]})
cases=[(broad,bb,BROAD,r'$S\sim U(0,60\,\mathrm{s})$','(a) Broad prior: long red remains plausible',r'$P(S<21.6\,\mathrm{s})=36\%$'),
       (late,lb,LATE,r'$S\sim U(26,30\,\mathrm{s})$','(b) Late prior: all switches require delay',r'$P(S<21.6\,\mathrm{s})=0$')]
for j,(pol,ub,col,plabel,title,note) in enumerate(cases):
    ax=axs[0,j]
    ax.plot(-kb['d'],kb['v']*3.6,color=GRAY,ls=(0,(6,3)),lw=2.4,label=r'Perfect knowledge: $S=28$ s')
    ax.plot(-ub['d'],ub['v']*3.6,color=col,lw=2.6,label='Imperfect knowledge: '+plabel)
    # exact switch markers
    for br,c in [(kb,GRAY),(ub,col)]:
        x=float(np.interp(28,br['t'],-br['d'])); v=float(np.interp(28,br['t'],br['v'])*3.6)
        ax.plot(x,v,'o',ms=6.4,color=c,mec='white',mew=1.1,zorder=8)
    ax.axvline(0,color=LIGHT,ls=':',lw=1.35)
    ax.set_title(title,fontsize=11.2,fontweight='semibold')
    ax.set_xlim(-305,105);ax.set_ylim(-2,53);ax.grid(alpha=.18);ax.spines[['top','right']].set_visible(False)
    ax.text(.03,.07,note,transform=ax.transAxes,fontsize=9.7,bbox=dict(fc='white',ec='none',alpha=.85,pad=2.5))
    if j==0:ax.set_ylabel('Speed (km/h)')
    ax.legend(loc='center left', bbox_to_anchor=(-295,15), bbox_transform=ax.transData, fontsize=8.1, framealpha=.9)
    at=axs[1,j]
    at.plot(kb['t'],-kb['d'],color=GRAY,ls=(0,(6,3)),lw=2.2)
    at.plot(ub['t'],-ub['d'],color=col,lw=2.4)
    at.axvline(28,color=LIGHT,ls=':',lw=1.25)
    at.axhline(0,color='#777',ls=':',lw=1.0)
    at.set_xlim(0,38);at.set_ylim(-305,105);at.grid(alpha=.18);at.spines[['top','right']].set_visible(False)
    at.set_xlabel('Time since first observing red (s)')
    if j==0:at.set_ylabel('Position relative to stop line (m)')
    if j==0:
        at.annotate('waiting at the line',xy=(26.4,0),xytext=(17.5,55),fontsize=8.5,
                    arrowprops=dict(arrowstyle='->',lw=1,color='#555'))
fig.suptitle('Same realized switch $S=28$ s; only prior information changes',fontsize=12.5,y=.995)
fig.text(.5,.005,'Dashed = perfect knowledge; solid = imperfect-information feedback; dots = green observed at 28 s',ha='center',fontsize=9.6)
fig.tight_layout(rect=[0,.025,1,.965])
fig.savefig(FIG/'fig1_information_comparison.png',dpi=250,bbox_inches='tight',facecolor='white');plt.close(fig)

# Fig. 2: conditional-red policies in position space, consistent colors.
fig,(ax,az)=plt.subplots(1,2,figsize=(12.2,4.9),gridspec_kw={'width_ratios':[2.05,1]})
for a in (ax,az):
    a.axvline(0,color=LIGHT,ls=':',lw=1.35);a.grid(alpha=.18);a.spines[['top','right']].set_visible(False);a.set_ylim(-1,53);a.set_xlabel('Position relative to stop line (m)')
ax.set_ylabel('Speed (km/h)')
for r,c,ls,label in [(known,GRAY,(0,(6,3)),r'Perfect knowledge: $S=28$ s'),(broad,BROAD,'-',r'Broad prior $U(0,60)$, conditional on red')]:
    x=-r['d'];v=r['v']*3.6
    ax.plot(x,v,color=c,ls=ls,lw=2.5,label=label);az.plot(x,v,color=c,ls=ls,lw=2.5)
    idx=np.where(r['b']>.95)[0]
    if len(idx):
        x0=float(x[idx[0]]);x1=float(x[min(idx[-1]+1,len(x)-1)])
        for a in (ax,az):a.axvspan(min(x0,x1),max(x0,x1),color=c,alpha=.09,zorder=0)
ax.set_xlim(-305,5);az.set_xlim(-25,2);az.set_xticks([-25,-20,-15,-10,-5,0])
ax.set_title('Full 300 m approach',fontsize=11);az.set_title('Last 25 m',fontsize=11)
ax.legend(loc='lower left',fontsize=8.7)
fig.suptitle('Conditional policy if the signal is still red',fontsize=12.5)
fig.text(.5,.015,'Light shading marks intervals in which full braking is active.',ha='center',fontsize=9.3)
fig.tight_layout(rect=[0,.04,1,.94])
fig.savefig(FIG/'fig2_red_conditional.png',dpi=250,bbox_inches='tight',facecolor='white');plt.close(fig)

# Figs. 3-4: narrow prior, 900 CZK/h, exact-cell-weight N=1100 data.
narrow=load('narrow900_N1100.npz',27.5,cost_key='J')
noinc=load('narrow900_noacc_N1100.npz',27.5,cost_key='J')
model900=Model(900)
active=(narrow['u']>.2)&(narrow['v'][:-1]<VMAX-.03)
inds=np.where(active)[0]; istart=inds[0]
t_acc=float(narrow['t'][istart]);x_acc=float(-narrow['d'][istart]);v_acc=float(narrow['v'][istart]*3.6)
branches=[]
for sw,c in [(27.03,'#4E82AC'),(27.25,GREEN),(27.47,'#8A5FA8')]:branches.append((sw,c,branch(narrow,sw,model900)))
fig,(ax,az)=plt.subplots(1,2,figsize=(12.5,5.4),gridspec_kw={'width_ratios':[1.85,1.05]})
for a in (ax,az):
    a.axvline(0,color=LIGHT,ls=':',lw=1.35);a.grid(alpha=.18);a.spines[['top','right']].set_visible(False);a.set_ylim(0,52);a.set_xlabel('Position relative to stop line (m)')
ax.set_ylabel('Speed (km/h)')
ax.plot(-noinc['d'],noinc['v']*3.6,color=COUNTER,ls=':',lw=1.7,label=r'Counterfactual: $u\leq u_{\rm cruise}$ while red')
az.plot(-noinc['d'],noinc['v']*3.6,color=COUNTER,ls=':',lw=1.7)
for sw,c,br in branches:
    for a in (ax,az):a.plot(-br['d'],br['v']*3.6,color=c,lw=1.7,alpha=.88)
    xsw=float(np.interp(sw,br['t'],-br['d']));vsw=float(np.interp(sw,br['t'],br['v'])*3.6)
    ax.scatter([xsw],[vsw],s=58,color=c,ec='white',lw=1,zorder=7,label=f'Green at {sw:.2f} s')
    az.scatter([xsw],[vsw],s=58,color=c,ec='white',lw=1,zorder=7)
for a in (ax,az):
    a.plot(-narrow['d'],narrow['v']*3.6,color=NARROW,lw=3.0,zorder=5,label=r'Narrow-prior policy while red' if a is ax else None)
    a.plot([x_acc],[v_acc],marker='D',color=NARROW,ms=6,zorder=8)
ax.set_xlim(-305,102);az.set_xlim(-17,3);az.set_xticks([-15,-10,-5,0])
ax.set_title('Whole trajectory',fontsize=11);az.set_title('Last 17 m',fontsize=11)
az.annotate('acceleration begins\nbefore green is possible',xy=(x_acc,v_acc),xytext=(-16,15),fontsize=8.8,
            arrowprops=dict(arrowstyle='->',lw=1,color='#444'),bbox=dict(fc='white',ec='none',alpha=.8,pad=1.4))
ax.legend(loc='lower left',fontsize=8.2,framealpha=.92)
fig.suptitle(r'Anticipatory acceleration with $S\sim U(27,27.5)$ s',fontsize=12.5)
fig.tight_layout(rect=[0,0,1,.95])
fig.savefig(FIG/'fig3_anticipatory_position.png',dpi=250,bbox_inches='tight',facecolor='white');plt.close(fig)

fig,axs=plt.subplots(2,1,figsize=(9.2,6.6),sharex=True,gridspec_kw={'height_ratios':[1.55,1]})
axs[0].plot(narrow['t'],narrow['v']*3.6,color=NARROW,lw=2.6,label='Narrow-prior policy while red')
axs[0].plot(noinc['t'],noinc['v']*3.6,color=COUNTER,ls=':',lw=1.7,label=r'Counterfactual: $u\leq u_{\rm cruise}$')
# extend step traces to the final time
xstep=narrow['t']; ustep=np.r_[narrow['u'],narrow['u'][-1]]; bstep=np.r_[narrow['b'],narrow['b'][-1]]
axs[1].step(xstep,ustep,where='post',color='#2F7B80',lw=2,label='Throttle $u$')
axs[1].step(xstep,bstep,where='post',color='#B0525C',lw=1.7,label='Brake $b$')
for a in axs:
    a.axvline(27.0,color='#7B8696',ls=(0,(6,3)),lw=1.3,label='Earliest possible green' if a is axs[0] else None)
    a.axvline(27.5,color='#7B8696',ls='-.',lw=1.3,label='Latest possible green' if a is axs[0] else None)
    a.grid(alpha=.18);a.spines[['top','right']].set_visible(False)
axs[0].plot([t_acc],[v_acc],marker='D',color=NARROW,ms=6,zorder=8)
axs[0].annotate('acceleration begins',xy=(t_acc,v_acc),xytext=(25.7,39.5),fontsize=8.8,
                arrowprops=dict(arrowstyle='->',lw=1,color='#444'))
axs[0].set_ylabel('Speed (km/h)');axs[0].set_ylim(20,48);axs[0].set_xlim(25.5,27.52);axs[0].legend(loc='upper left',fontsize=8.5)
axs[1].set_ylabel('Control (0-1)');axs[1].set_xlabel('Time since first observing red (s)');axs[1].set_ylim(-.06,1.06);axs[1].legend(loc='upper left',fontsize=8.5)
fig.suptitle('Acceleration begins before green is possible',fontsize=12.5)
fig.tight_layout(rect=[0,0,1,.96])
fig.savefig(FIG/'fig4_anticipatory_controls.png',dpi=250,bbox_inches='tight',facecolor='white');plt.close(fig)
print('Wrote figures to',FIG)

# Render branch-aware Fig. 1 last so it does not get overwritten by the older two-policy plot.
import runpy
runpy.run_path(str(ROOT/'make_fig1_unambiguous.py'),run_name='__main__')
