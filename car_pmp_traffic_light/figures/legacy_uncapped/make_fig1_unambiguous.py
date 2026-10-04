from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'numerics'))
from model import Model
from solve_policies import branch

DATA=ROOT/'data'
OUT=ROOT/'figures'
GRAY='#303641'; BROWN='#BF6A36'; BROAD=BROWN; LATE=BROWN; LIGHT='#A23B3B'; GREEN='#148267'

def load(name,T):
    a=np.load(DATA/name)
    r={k:a[k] for k in a.files}; r['T']=T
    return r

base=Model(300)
known=load('k_28.npz',28)
broad=load('u_0_60.npz',60)
late=load('u_26_30.npz',30)
kb=branch(known,28,base)
bb=branch(broad,28,base)
late_branches={S:branch(late,S,base) for S in [26.,27.,28.,29.,30.]}
lb28=late_branches[28.]

fig,axs=plt.subplots(2,2,figsize=(12.2,7.2),gridspec_kw={'height_ratios':[1.25,1.0]})

# --- (a) broad prior: common still-red branch plus the realized S=28 continuation ---
ax=axs[0,0]
ax.plot(-kb['d'],kb['v']*3.6,color=GRAY,ls=(0,(6,3)),lw=2.4,zorder=3,
        label=r'Perfect knowledge: $S=28$ s')
# Thick brown is strictly the conditional policy if red continues to be observed.
ax.plot(-broad['d'],broad['v']*3.6,color=BROAD,lw=3.0,zorder=4,
        label=r'Imperfect: $U(0,60)$, still red')
# The realized S=28 continuation leaves the red branch only after green is observed.
mask=bb['t']>=28-1e-9
ax.plot(-bb['d'][mask],bb['v'][mask]*3.6,color=GREEN,lw=2.5,alpha=1.0,zorder=6,
        label=r'After observing green at $28$ s')
xs=float(np.interp(28,bb['t'],-bb['d'])); vs=float(np.interp(28,bb['t'],bb['v'])*3.6)
ax.plot(xs,vs,'o',ms=6.2,color=GREEN,mec='white',mew=1.0,zorder=8)
# Mark the known-S=28 switch on the dashed reference as well.
xk=float(np.interp(28,kb['t'],-kb['d'])); vk=float(np.interp(28,kb['t'],kb['v'])*3.6)
ax.plot(xk,vk,'o',ms=5.6,color=GRAY,mec='white',mew=1.0,zorder=8)
ax.axvline(0,color=LIGHT,ls=':',lw=1.35)
ax.set_title('(a) Broad prior: long red remains plausible',fontsize=11.2,fontweight='semibold')
ax.set_xlim(-305,105);ax.set_ylim(-2,53);ax.grid(alpha=.18);ax.spines[['top','right']].set_visible(False)
ax.text(.03,.07,r'$P(S<21.6\,\mathrm{s})=36\%$',transform=ax.transAxes,fontsize=9.7,bbox=dict(fc='white',ec='none',alpha=.85,pad=2.5))
ax.set_ylabel('Speed (km/h)')
ax.legend(loc='center left', bbox_to_anchor=(-295,15), bbox_transform=ax.transData, fontsize=7.9, framealpha=.9)

at=axs[1,0]
at.plot(kb['t'],-kb['d'],color=GRAY,ls=(0,(6,3)),lw=2.2,zorder=3)
at.plot(broad['t'],-broad['d'],color=BROAD,lw=3.0,zorder=4)
at.plot(bb['t'][mask],-bb['d'][mask],color=GREEN,lw=2.5,alpha=1.0,zorder=6)
at.plot(28,xs,'o',ms=6.2,color=GREEN,mec='white',mew=1.0,zorder=8)
at.axvline(28,color=LIGHT,ls=':',lw=1.25)
at.axhline(0,color='#777',ls=':',lw=1.0)
at.set_xlim(0,38);at.set_ylim(-305,105);at.grid(alpha=.18);at.spines[['top','right']].set_visible(False)
at.set_xlabel('Time since first observing red (s)');at.set_ylabel('Position relative to stop line (m)')
at.annotate('still red: waiting at the line',xy=(27.0,0),xytext=(16.7,55),fontsize=8.5,
            arrowprops=dict(arrowstyle='->',lw=1,color='#555'))

# --- (b) late prior, branch-aware ---
ax=axs[0,1]
# perfect-information reference
ax.plot(-kb['d'],kb['v']*3.6,color=GRAY,ls=(0,(6,3)),lw=2.4,zorder=3)
# common red branch only up to latest green bound; this is the feedback policy before observation
ax.plot(-late['d'],late['v']*3.6,color=LATE,lw=3.0,zorder=4)
# All post-switch branches are GREEN and thin. Brown is strictly still-red.
# Importantly, S=30 remains on the brown branch until x=0 / t=30.
for S,br in late_branches.items():
    mask_b=br['t']>=S-1e-9
    is28=abs(S-28.0)<1e-9
    ax.plot(-br['d'][mask_b], br['v'][mask_b]*3.6, color=GREEN,
            lw=2.5 if is28 else 1.0, alpha=1.0 if is28 else 0.45,
            zorder=6 if is28 else 2)
    xs_b=float(np.interp(S,br['t'],-br['d'])); vs_b=float(np.interp(S,br['t'],br['v'])*3.6)
    ax.plot(xs_b,vs_b,'o',ms=6.2 if is28 else 4.1,color=GREEN,
            mec='white',mew=1.0 if is28 else 0.8,
            alpha=1.0 if is28 else 0.72,zorder=8)
ax.axvline(0,color=LIGHT,ls=':',lw=1.35)
ax.set_title('(b) Late prior: common red branch and green futures',fontsize=11.2,fontweight='semibold')
ax.set_xlim(-305,105);ax.set_ylim(-2,53);ax.grid(alpha=.18);ax.spines[['top','right']].set_visible(False)
ax.text(.03,.07,r'$P(S<21.6\,\mathrm{s})=0$',transform=ax.transAxes,fontsize=9.7,bbox=dict(fc='white',ec='none',alpha=.85,pad=2.5))
# custom legend, intentionally at ~15 km/h left
handles=[
    Line2D([0],[0],color=GRAY,ls=(0,(6,3)),lw=2.4,label=r'Perfect knowledge: $S=28$ s'),
    Line2D([0],[0],color=LATE,lw=3.0,label=r'Imperfect: $U(26,30)$, still red'),
    Line2D([0],[0],color=GREEN,lw=2.5,marker='o',mec='white',label=r'After green ($28$ s emphasized)')]
ax.legend(handles=handles,loc='center left', bbox_to_anchor=(-295,15), bbox_transform=ax.transData, fontsize=8.0, framealpha=.9)
# Small inset: makes the timing distinction explicit even in a two-column PDF.
iax=ax.inset_axes([.58,.045,.39,.35])
iax.plot(-kb['d'],kb['v']*3.6,color=GRAY,ls=(0,(6,3)),lw=1.5)
iax.plot(-late['d'],late['v']*3.6,color=LATE,lw=2.6)
for S,br in late_branches.items():
    mask_b=br['t']>=S-1e-9
    is28=abs(S-28.0)<1e-9
    iax.plot(-br['d'][mask_b],br['v'][mask_b]*3.6,color=GREEN,
             lw=1.8 if is28 else .8,alpha=1.0 if is28 else .45)
    xs_b=float(np.interp(S,br['t'],-br['d']));vs_b=float(np.interp(S,br['t'],br['v'])*3.6)
    iax.plot(xs_b,vs_b,'o',ms=3.8 if is28 else 2.6,color=GREEN,
             mec='white',mew=.6,alpha=1.0 if is28 else .7)
iax.axvline(0,color=LIGHT,ls=':',lw=1)
iax.set_xlim(-39,13);iax.set_ylim(26,52)
iax.set_title('Near the stop line',fontsize=7.6,pad=2)
iax.tick_params(labelsize=7,length=2)
iax.grid(alpha=.15)


at=axs[1,1]
at.plot(kb['t'],-kb['d'],color=GRAY,ls=(0,(6,3)),lw=2.2,zorder=3)
at.plot(late['t'],-late['d'],color=LATE,lw=3.0,zorder=4)
for S,br in late_branches.items():
    mask_b=br['t']>=S-1e-9
    is28=abs(S-28.0)<1e-9
    at.plot(br['t'][mask_b],-br['d'][mask_b],color=GREEN,
            lw=2.5 if is28 else 1.0,alpha=1.0 if is28 else .45,
            zorder=6 if is28 else 2)
    xs_b=float(np.interp(S,br['t'],-br['d']))
    at.plot(S,xs_b,'o',ms=6.2 if is28 else 4.1,
            color=GREEN,mec='white',mew=1.0 if is28 else .8,zorder=8,
            alpha=1.0 if is28 else .72)
    offsets={26:(3,6),27:(3,-13),28:(3,6),29:(3,-13),30:(3,6)}
    dx,dy=offsets[S]
    at.annotate(f'{S:g}',(S,xs_b),xytext=(dx,dy),
                textcoords='offset points',fontsize=7.8 if is28 else 7.3,
                fontweight='semibold' if is28 else 'normal',color=GREEN,
                alpha=1.0 if is28 else .75)
at.axvline(28,color=LIGHT,ls=':',lw=1.0,alpha=.7)
at.axhline(0,color='#777',ls=':',lw=1.0)
at.set_xlim(0,38);at.set_ylim(-305,105);at.grid(alpha=.18);at.spines[['top','right']].set_visible(False)
at.set_xlabel('Time since first observing red (s)')
at.text(.52,.09,'Numbers = observed green time (s)',transform=at.transAxes,fontsize=8.3,
        bbox=dict(fc='white',ec='none',alpha=.82,pad=2.2))

fig.suptitle('Perfect knowledge versus feedback: the common red branch and possible green continuations',fontsize=12.3,y=.995)
fig.text(.5,.005,'Dashed = perfect knowledge; thick brown = still red; green = after green observed; the 28 s green branch is emphasized',ha='center',fontsize=9.5)
fig.tight_layout(rect=[0,.025,1,.965])

png=OUT/'fig1_information_comparison.png'
pdf=OUT/'fig1_information_comparison.pdf'
fig.savefig(png,dpi=300,bbox_inches='tight',facecolor='white')
fig.savefig(pdf,bbox_inches='tight',facecolor='white')
plt.close(fig)
print(png)
print(pdf)
