"""Session 5B publication plots: capped acceleration, 4 m/s2 braking.

The article and its historical figure paths remain unmodified until Session 6.
Figures made by this script are `fig{1..4}_*_capped` with *distinct* names.

Brown is reserved solely for a common conditional trajectory while red is
still observed. Green (and other branch colors) begin at an observed switch;
dashed charcoal denotes the perfect-information deterministic reference.
"""
from pathlib import Path
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'numerics'))
from model import Model, D0, VMAX, published_accel_cap
from solve_policies import branch

FIG = ROOT/'figures'
DATA = ROOT/'data'
FIG.mkdir(exist_ok=True)
BROWN = '#A85C30'  # same semantic brown; stronger contrast with green
GREEN = '#147C62'
GREEN_FAINT = '#61A791'
GREEN_BLUE = '#497DA6'  # permissible branch distinction
GREEN_PURPLE = '#8865A1'
DARK = '#344052'
SLATE = '#757F90'
RED = '#B04747'
TEAL = '#267E8A'
GRID = '#CBD2D9'
PURPLE = '#7560A1'

plt.rcParams.update({
    'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':10.5,
    'axes.labelsize':10,'xtick.labelsize':9,'ytick.labelsize':9,
    'legend.fontsize':8.6,'pdf.fonttype':42,'ps.fonttype':42,
    'figure.facecolor':'white','savefig.facecolor':'white',
})

M300=Model(300,accel_cap=published_accel_cap(),brake_a=4,rest_mode='balance')
M900=Model(900,accel_cap=published_accel_cap(),brake_a=4,rest_mode='balance')


def read(path, required_model=False):
    with np.load(path,allow_pickle=False) as z:
        p={k:z[k].copy() for k in z.files}
    if required_model:
        assert str(p['model_id'])=='capped_b4_balance',path
        assert abs(float(p['accel_cap'])-M300.accel_cap)<1e-9
        assert abs(float(p['brake_a'])-4.0)<1e-9
    assert np.all(np.diff(p['t'])>0) and np.all(np.isfinite(p['v']))
    assert len(p['u'])==len(p['t'])-1==len(p['b'])
    return p


def format_axis(ax, xlab=None, ylab=None, xlim=None, ylim=None):
    ax.spines[['top','right']].set_visible(False)
    ax.grid(color=GRID,lw=.6,alpha=.65,zorder=0)
    if xlab:ax.set_xlabel(xlab)
    if ylab:ax.set_ylabel(ylab)
    if xlim:ax.set_xlim(*xlim)
    if ylim:ax.set_ylim(*ylim)


def save(fig,name):
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(FIG/(name+'.png'),dpi=210,bbox_inches='tight')
    plt.close(fig)


def only_green(ax,policy,S,model,**kwargs):
    q=branch(policy,S,model)
    after=q['t']>=S-1e-9
    ax.plot(-q['d'][after],q['v'][after]*3.6,**kwargs)
    return q


def green_time(ax,policy,S,model,which='pos',**kwargs):
    q=branch(policy,S,model)
    after=q['t']>=S-1e-9
    yy = (-q['d'] if which=='pos' else q['v']*3.6)
    ax.plot(q['t'][after],yy[after],**kwargs)
    return q


# Figure 1: deterministic known-time mechanism, exact reconstructed N=480.
known=read(ROOT/'diagnostics'/'capped_known_sweep'/'capped_b4_028.00_N480.npz',True)
assert abs(float(known['cost'])-3.6755043875231594)<1e-7
S=28.;t=known['t'];u=known['u'];b=known['b'];speed=known['v']*3.6
nup=np.flatnonzero((u>.2)&(known['v'][:-1]<VMAX-.03))
onset=float(t[nup[0]])
fig=plt.figure(figsize=(12.1,5.7),constrained_layout=False)
gs=fig.add_gridspec(2,2,width_ratios=[1.75,1],height_ratios=[1.1,.7],left=.07,right=.985,top=.9,bottom=.13,wspace=.22,hspace=.3)
a=fig.add_subplot(gs[0,0]);ac=fig.add_subplot(gs[1,0],sharex=a);ap=fig.add_subplot(gs[:,1])
a.plot(t,speed,color=DARK,lw=2.6,label=r'Known green at $S=28$ s')
a.axvline(28,color=RED,ls=':',lw=1.5)
a.scatter([28],[speed[-1]],color=DARK,edgecolors='white',zorder=6,s=48)
a.annotate('Speed recovery begins',xy=(onset,np.interp(onset,t,speed)),xytext=(16,45),
           fontsize=9,arrowprops={'arrowstyle':'->','color':DARK,'lw':1.1})
a.text(.04,.06,'Cruise-to-line: 21.6 s\nFree-coast-to-line: 25.3 s\nRequired switch: 28 s',transform=a.transAxes,
       fontsize=9,va='bottom',bbox={'fc':'white','ec':'#D1D5DB','boxstyle':'round,pad=.35'})
format_axis(a,ylab='Speed (km/h)',xlim=(0,28.4),ylim=(27,53))
a.set_title('(a) Speed: delay first, then recover',loc='left',fontweight='semibold')
ac.step(t,np.r_[u,u[-1]],where='post',color=TEAL,lw=1.8,label='Throttle')
ac.step(t,np.r_[b,b[-1]],where='post',color=RED,lw=1.8,label='Brake')
ac.axvline(28,color=RED,ls=':',lw=1.25)
ac.annotate('Initial braking',xy=(.14,1),xytext=(3.3,.79),fontsize=8.7,
            arrowprops={'arrowstyle':'->','color':RED,'lw':.9})
format_axis(ac,xlab='Time since first red observation (s)',ylab='Control (0–1)',xlim=(0,28.4),ylim=(-.07,1.12))
ac.legend(loc='center left',bbox_to_anchor=(.32,.8),frameon=False,ncol=2)
ac.set_title('(b) Brake–coast–accelerate control',loc='left',fontweight='semibold')
ap.plot(-known['d'],speed,color=DARK,lw=2.5)
ap.axvline(0,color=RED,ls=':',lw=1.6)
ap.plot([0],[speed[-1]],'o',color=DARK,ms=5)
ap.scatter([-np.interp(onset,t,known['d'])],[np.interp(onset,t,speed)],s=51,marker='D',color=TEAL,zorder=5)
format_axis(ap,xlab='Position relative to stop line (m)',ylab='Speed (km/h)',xlim=(-65,4),ylim=(27,52))
ap.set_title('(c) Recovery near the stop line',loc='left',fontweight='semibold')
fig.suptitle('Known switch at 28 s: braking supplies delay, acceleration preserves speed at green',fontsize=13,fontweight='semibold',y=.985)
fig.text(.50,.035,'Capped acceleration (4.137 m/s²), maximum brake parameter 4 m/s²; red dotted lines mark the stop line or switch.',ha='center',fontsize=9,color='#454E5B')
save(fig,'fig1_known28_capped')

# Figure 2: partial information with common-red trajectory and contingent green branches.
broad=read(DATA/'u_0_60_capped.npz',True)
late=read(DATA/'u_26_30_capped.npz',True)
fig,axs=plt.subplots(2,2,figsize=(12.4,7.2),gridspec_kw={'height_ratios':[1.17,1.]})
for j,(policy,prior,subtitle) in enumerate([
  (broad,r'Broad prior  $S\sim U(0,60)$','Early green possible; late red may continue'),
  (late,r'Late prior  $S\sim U(26,30)$','Delay required; approach remains rolling'),
]):
    aa,pt=axs[0,j],axs[1,j]
    br_known=branch(known,28,M300)
    aa.plot(-br_known['d'],br_known['v']*3.6,color=DARK,lw=2,ls=(0,(5,3)))
    pt.plot(br_known['t'],-br_known['d'],color=DARK,lw=2,ls=(0,(5,3)))
    aa.plot(-policy['d'],policy['v']*3.6,color=BROWN,lw=2.75,zorder=4)
    pt.plot(policy['t'],-policy['d'],color=BROWN,lw=2.75,zorder=4)
    for s in ((28.,) if j==0 else (26.,28.,30.)):
        focus=s==28.
        color=GREEN
        br=only_green(aa,policy,s,M300,color=color,lw=2.6 if focus else 1.05,alpha=1. if focus else .42,zorder=6 if focus else 2)
        green_time(pt,policy,s,M300,color=color,lw=2.6 if focus else 1.05,alpha=1. if focus else .42,zorder=6 if focus else 2)
        sx=float(np.interp(s,br['t'],-br['d']))
        sy=float(np.interp(s,br['t'],br['v'])*3.6)
        aa.scatter([sx],[sy],color=color,s=54 if focus else 27,edgecolors='white',linewidths=.8,zorder=8)
        pt.scatter([s],[sx],color=color,s=54 if focus else 27,edgecolors='white',linewidths=.8,zorder=8)
        if j==1:
            pt.annotate(f'{s:g}',xy=(s,sx),xytext=(3,5 if s!=28 else -14),textcoords='offset points',fontsize=8.5,color=GREEN)
    aa.axvline(0,color=RED,ls=':',lw=1.2)
    pt.axhline(0,color=RED,ls=':',lw=1.2)
    pt.axvline(28,color='#C36D64',ls=':',lw=1.0)
    format_axis(aa,xlab='Position relative to stop line (m)',ylab='Speed (km/h)' if j==0 else None,xlim=(-305,104),ylim=(-2,53))
    format_axis(pt,xlab='Time since first red observation (s)',ylab='Position relative to line (m)' if j==0 else None,xlim=(0,37),ylim=(-305,105))
    aa.set_title('('+('a' if j==0 else 'b')+') '+prior,loc='left',fontweight='semibold')
    pt.set_title(subtitle,loc='left',fontsize=10)
    if j==0:
        pt.annotate('Stop and wait if red persists',xy=(28,0),xytext=(16,62),fontsize=9,
                    arrowprops={'arrowstyle':'->','lw':1.1,'color':BROWN})
legend=[Line2D([0],[0],color=DARK,ls=(0,(5,3)),lw=2,label='Known $S=28$ s'),
        Line2D([0],[0],color=BROWN,lw=2.8,label='Common red-conditional policy'),
        Line2D([0],[0],color=GREEN,lw=2.7,marker='o',label='Green observed; $S=28$ emphasized')]
fig.legend(handles=legend,loc='lower center',bbox_to_anchor=(.5,.024),ncol=3,frameon=False,fontsize=9)
fig.suptitle('Same realized green at 28 s: information changes when delay is created',fontsize=13,fontweight='semibold',y=.99)
fig.tight_layout(rect=[.0,.065,1,.95],h_pad=1.8,w_pad=1.6)
save(fig,'fig2_information_capped')

# Figure 3: narrowly uncertain switch; causally valid before-27 comparison.
unrestricted=read(DATA/'s4_capped900_unrestricted_N1100.npz')
pre27=read(DATA/'s4_capped900_before27_N1100.npz')
all_red=read(DATA/'s4_capped900_all_red_N1100.npz')
assert np.allclose(unrestricted['t'],pre27['t'])
assert abs(float(unrestricted['cost'])-9.40022255399851)<1e-6
adv=float(pre27['cost'])-float(unrestricted['cost'])
assert abs(adv-.010246485626287)<1e-6
# onset = first speed-increasing throttle rather than u > arbitrary threshold
uc=M900.cruise_throttle(unrestricted['v'][:-1])
mask=(unrestricted['u']>uc+.03)&(unrestricted['v'][:-1]<VMAX-.03)
start=int(np.flatnonzero(mask)[0]);ta=float(unrestricted['t'][start]);va=float(unrestricted['v'][start]*3.6);da=float(-unrestricted['d'][start])
assert ta<27.0
fig=plt.figure(figsize=(12.2,6.0))
gs=fig.add_gridspec(2,2,width_ratios=[1.65,1],height_ratios=[1.1,.72],left=.075,right=.985,top=.88,bottom=.13,hspace=.29,wspace=.24)
aspeed=fig.add_subplot(gs[0,0]);actrl=fig.add_subplot(gs[1,0],sharex=aspeed);apos=fig.add_subplot(gs[:,1])
aspeed.plot(unrestricted['t'],unrestricted['v']*3.6,color=BROWN,lw=2.8,label='Optimal candidate: still red')
aspeed.plot(pre27['t'],pre27['v']*3.6,color=SLATE,lw=2,ls=(0,(4,2)),label='No acceleration before 27 s')
aspeed.scatter([ta],[va],color=BROWN,marker='D',s=52,zorder=8)
aspeed.annotate(f'Acceleration starts near {ta:.1f} s',xy=(ta,va),xytext=(26.02,40.3),fontsize=9,
                arrowprops={'arrowstyle':'->','color':BROWN,'lw':1.15})
actrl.step(unrestricted['t'],np.r_[unrestricted['u'],unrestricted['u'][-1]],where='post',color=TEAL,lw=2,label='Throttle (unrestricted)')
actrl.step(pre27['t'],np.r_[pre27['u'],pre27['u'][-1]],where='post',color=SLATE,ls='--',lw=1.6,label='Throttle (restricted)')
actrl.step(unrestricted['t'],np.r_[unrestricted['b'],unrestricted['b'][-1]],where='post',color=RED,lw=1.1,alpha=.8,label='Brake')
for ax in (aspeed,actrl):
    ax.axvspan(25.8,27,facecolor=BROWN,alpha=.07,zorder=-1)
    ax.axvline(27,color=DARK,ls=(0,(5,3)),lw=1.15)
    ax.axvline(27.5,color=DARK,ls=':',lw=1.15)
format_axis(aspeed,ylab='Speed (km/h)',xlim=(25.8,27.52),ylim=(31,47))
format_axis(actrl,xlab='Time since first red observation (s)',ylab='Control (0–1)',xlim=(25.8,27.52),ylim=(-.04,1.07))
aspeed.set_title('(a) Speed before and during uncertain green window',loc='left',fontweight='semibold')
actrl.set_title('(b) Control: before-27 s restriction only',loc='left',fontweight='semibold')
aspeed.legend(loc='upper left',frameon=False,fontsize=8.3)
actrl.legend(loc='upper left',frameon=False,fontsize=7.9,ncol=1)
apos.plot(-unrestricted['d'],unrestricted['v']*3.6,color=BROWN,lw=2.8,label='Red-conditional trajectory',zorder=5)
apos.plot(-pre27['d'],pre27['v']*3.6,color=SLATE,lw=1.8,ls=(0,(4,2)),label='No acceleration before 27 s')
for s,color in [(27.03,GREEN_BLUE),(27.25,GREEN),(27.47,GREEN_PURPLE)]:
    br=only_green(apos,unrestricted,s,M900,color=color,lw=1.9,label=f'Green at {s:.2f} s')
    apos.scatter([np.interp(s,br['t'],-br['d'])],[np.interp(s,br['t'],br['v'])*3.6],s=44,color=color,edgecolors='white',linewidths=.8,zorder=8)
apos.plot(da,va,marker='D',ms=6,color=BROWN,zorder=9)
apos.axvline(0,color=RED,ls=':',lw=1.1)
format_axis(apos,xlab='Position relative to stop line (m)',ylab='Speed (km/h)',xlim=(-23,42),ylim=(30,52))
apos.set_title('(c) Multiple possible green continuations',loc='left',fontweight='semibold')
apos.legend(loc='lower right',frameon=True,framealpha=.91,fontsize=7.9)
fig.suptitle(r'Acceleration before green is possible, $S\sim U(27,27.5)$ s; time value 900 CZK/h',fontsize=12.5,fontweight='semibold',y=.98)
fig.text(.51,.041,fr'Expected cost: unrestricted {float(unrestricted["cost"]):.6f} CZK; no acceleration before 27 s {float(pre27["cost"]):.6f} CZK; benefit {adv:.6f} CZK.',
         ha='center',fontsize=9.1)
save(fig,'fig3_anticipation_capped')

# Figure 4: split axes because late prior is only four seconds wide.
with np.load(DATA/'voi_capped_regret.npz',allow_pickle=False) as z:
    assert str(z['model_id'])=='capped_b4_balance'
    rg={k:z[k].copy() for k in z.files}
with open(ROOT/'diagnostics'/'session5a_voi.json',encoding='utf-8') as f:
    numbers=json.load(f)
fig,axes=plt.subplots(1,2,figsize=(12.1,4.9),gridspec_kw={'width_ratios':[1.35,1]})
for j,(ax,s,r,title,key,color) in enumerate([
    (axes[0],rg['S'],rg['regret_broad'],r'(a) Broad prior: $U(0,60)$','U_0_60',DARK),
    (axes[1],rg['late_S'],rg['late_regret'],r'(b) Late prior: $U(26,30)$','U_26_30',PURPLE),
]):
    ax.plot(s,r,color=color,lw=2.45)
    ax.fill_between(s,0,r,color=color,alpha=.13)
    ax.axhline(0,color='#454545',lw=.8)
    xlims=(0,60) if j==0 else (26,30)
    format_axis(ax,xlab='Realized switch time S (s)',ylab='Regret (CZK)' if j==0 else None,
                xlim=xlims,ylim=(-.012,.55 if j==0 else .40))
    ax.set_title(title,loc='left',fontweight='semibold')
    vo=numbers[key]['value_information_CZK']
    ax.text(.045,.93,f'Average VOI = {vo:.6f} CZK',transform=ax.transAxes,
            fontsize=10,fontweight='semibold',va='top',bbox={'fc':'white','ec':'#D1D5DB','boxstyle':'round,pad=.4'})
fig.suptitle('Value of switch-time information: realized regret and prior-weighted average',fontsize=12.8,fontweight='semibold',y=.99)
fig.text(.5,.035,r'$R(s)=J_{\mathrm{feedback}}(s)-J_{\mathrm{known}}^*(s)$; shaded area / prior-support width gives average regret (VOI).',ha='center',fontsize=10)
fig.tight_layout(rect=[0,.065,1,.92],w_pad=2.)
save(fig,'fig4_regret_capped')

manifest={
 'model_id':'capped_b4_balance','accel_cap_m_s2':M300.accel_cap,'brake_parameter_m_s2':4.0,
 'rest_mode':'balance','deterministic_known_S_28_cost_CZK':float(known['cost']),
 'deterministic_first_speed_increasing_throttle_s':onset,
 'broad_feedback_expected_cost_CZK':float(broad['cost']),
 'late_feedback_expected_cost_CZK':float(late['cost']),
 'anticipation_pre27_acceleration_onset_s':ta,
 'anticipation_cost_unrestricted_CZK':float(unrestricted['cost']),
 'anticipation_cost_no_pre27_acceleration_CZK':float(pre27['cost']),
 'anticipation_pre27_benefit_CZK':adv,
 'broad_VOI_CZK':numbers['U_0_60']['value_information_CZK'],
 'late_VOI_CZK':numbers['U_26_30']['value_information_CZK'],
 'files':['fig1_known28_capped','fig2_information_capped','fig3_anticipation_capped','fig4_regret_capped'],
 'numerical_status':'best-found feasible local candidates, not global certificates',
 'paper_status':'article.tex and article.pdf are not updated until Session 6',
}
(FIG/'session5b_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Saved four publication figure pairs (PDF+PNG) and figure manifest; capped model only.')
print(json.dumps({k:v for k,v in manifest.items() if k.endswith('_CZK') or k.endswith('_s')},indent=2))
