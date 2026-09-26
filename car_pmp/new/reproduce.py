"""Reproduce the article's fit, numerical results and figures (Python 3).
Dependencies: numpy, scipy, matplotlib. Run from any working directory.
All dynamics use metres and seconds; reported costs use CZK.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
(ROOT/'figures').mkdir(exist_ok=True)
Vdata=np.array([9.,80.,130.]); FCdata=np.array([10.,5.,6.5])
vdata=Vdata/3.6
Q0=0.6
c1,c2,c3=np.linalg.solve(np.column_stack([vdata,vdata**2,vdata**3]),FCdata*Vdata/100-Q0)
def Q(v): return Q0+c1*np.asarray(v)+c2*np.asarray(v)**2+c3*np.asarray(v)**3
def Qp(v): return c1+2*c2*v+3*c3*v*v
def FC(V): return 100*Q(np.asarray(V)/3.6)/V

m=1500.; Pe=110000.; eta_d=.9; Pw=eta_d*Pe
Cd=.302; area=2.2; rho=1.225; Crr=.012; grav=9.81; v0=1.; a=7.
A=Pw/m; B=rho*Cd*area/(2*m); r0=Crr*grav
ct=300/3600; cf=32.; qi=Q0/3600
vmin=60/3.6; vmax=160/3.6
L=1000.; X0=-2000.; X1=3000.
limits={'out':(130/3.6,160/3.6),'in':(80/3.6,110/3.6)}
phi_low=65/3600; phi_high=390/3600

def R(v): return r0+B*np.asarray(v)**2
def h(v): return A/(np.asarray(v)+v0)
def F(v): return h(v)-R(v)
def ucr(v): return R(v)/h(v)
def Gamma(v): return (Q(v)/3600-qi)*h(v)/R(v)
def phi(v,reg):
 l,t=limits[reg]
 z=np.asarray(v)
 return np.where(z<=l+1e-12,0.,np.where(z<=t+1e-12,phi_low,phi_high))
def g(v,reg): return (ct+cf*Q(v)/3600+phi(v,reg))/v

def cruise_candidates(reg):
 bounds=sorted(set([vmin,vmax]+[x for x in limits[reg] if vmin<x<vmax]))
 candidates=list(bounds)
 for low,up in zip(bounds[:-1],bounds[1:]):
  p=float(phi((low+up)/2,reg))
  s=lambda v: cf*(v*Qp(v)-Q(v))/3600-ct-p
  if s(low)*s(up)<0: candidates.append(brentq(s,low,up))
 candidates=[v for v in candidates if ucr(v)<=1+1e-12]
 return sorted([(v*3.6,float(g(v,reg)*1000)) for v in candidates],key=lambda x:x[1])
vo=cruise_candidates('out')[0][0]/3.6
vi=cruise_candidates('in')[0][0]/3.6
Go=float(g(vo,'out')); Gi=float(g(vi,'in'))
thresholds=sorted(set(x for ls in limits.values() for x in ls))
def integ(fun,lo,hi):
 if abs(hi-lo)<1e-12:return 0.
 return quad(fun,lo,hi,points=[x for x in thresholds if lo<x<hi],epsabs=1e-10)[0]

def parts(vb,vc,brake):
 rd=lambda v:R(v)+a*brake
 dout=integ(lambda v:v/rd(v),vb,vo); din=integ(lambda v:v/rd(v),vi,vb)
 ain=integ(lambda v:v/F(v),vi,vc); aout=integ(lambda v:v/F(v),vc,vo)
 if din+ain>L+1e-8:return None
 jd=integ(lambda v:(ct+cf*qi+phi(v,'out'))/rd(v),vb,vo)+integ(lambda v:(ct+cf*qi+phi(v,'in'))/rd(v),vi,vb)
 ja=integ(lambda v:(ct+cf*(qi+Gamma(v))+phi(v,'in'))/F(v),vi,vc)+integ(lambda v:(ct+cf*(qi+Gamma(v))+phi(v,'out'))/F(v),vc,vo)
 cost=Go*((-X0-dout)+(X1-L-aout))+jd+ja+Gi*(L-din-ain)
 return dict(cost=float(cost),Ddec=dout+din,Dacc=ain+aout,x12=-dout,x23=din,x34=L-ain,x45=L+aout)

# For this example N(v)>0 throughout (vi,vo): both variables minimize at vi.
# This certifies the boundary-speed minimum for either deceleration variant.
Nmid=(Go-Gi)*(vi+vo)/2+float(phi((vi+vo)/2,'in')-phi((vi+vo)/2,'out'))
assert Nmid>0 and (Go-Gi)*vo+phi_high>0
assert max(ucr(np.array([vi,vo])))<1 and min(F(np.array([vi,vo])))>0
rows={}
for label,brake in [('coast',0),('brake',1)]:
 p=parts(vi,vi,brake)
 rd=lambda v:R(v)+a*brake
 tc=integ(lambda v:1/rd(v),vi,vo)
 ta=integ(lambda v:1/F(v),vi,vo)
 outside=(X1-X0-L-p['Ddec']-p['Dacc'])
 t=outside/vo+tc+L/vi+ta
 fuel=outside/vo*Q(vo)/3600+qi*tc+L/vi*Q(vi)/3600+integ(lambda v:(qi+Gamma(v))/F(v),vi,vo)
 fine=float(phi(vi,'in'))*L/vi
 assert abs(p['cost']-(ct*t+cf*fuel+fine))<1e-8
 K=lambda v:cf*Gamma(v)/h(v)
 identity=Go*outside+Gi*L+integ(lambda v:g(v,'out')*v/rd(v),vi,vo)+integ(lambda v:g(v,'out')*v/F(v),vi,vo)+integ(lambda v:K(v)*a*brake/rd(v),vi,vo)
 assert abs(p['cost']-identity)<1e-8
 p.update(time=t,fuel=fuel,fines=fine,dec_time=tc,acc_time=ta)
 rows[label]=p
baseline=dict(time=(X1-X0)/vo,fuel=(X1-X0)/vo*Q(vo)/3600)
baseline['fines']=float(phi(vo,'in'))*L/vo
baseline['cost']=ct*baseline['time']+cf*baseline['fuel']+baseline['fines']
rows['constant_130']=baseline
results={'fuel_coefficients':[c1,c2,c3],'A':A,'B':B,'r0':r0,'cruise_candidates':{r:cruise_candidates(r) for r in ['out','in']},'rows':rows,'coast_130_120_seconds':integ(lambda v:1/R(v),120/3.6,130/3.6),'fuel_at_data':FC(Vdata).tolist(),'steady_130_power_kw':m*R(vo)*vo/1000,'steady_130_efficiency':m*R(vo)*vo/(Q(vo)*35.8e6/3600),'full_throttle_130_Lh':(qi+Gamma(vo))*3600,'model_terminal_speed_kmh':brentq(F,vo,100)*3.6}
(ROOT/'results.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2,'savefig.bbox':'tight'})
fig,ax=plt.subplots(figsize=(7.1,3.5))
V=np.linspace(6,160,600)
ax.plot(V,FC(V),color='#126782',lw=2,label='Empirical interpolation / extrapolation')
ax.errorbar(Vdata,FCdata,yerr=.5,fmt='o',color='#e07a24',capsize=4,label='Owner measurements; assumed ±0.5 tolerance')
ax.axvspan(130,160,color='#888888',alpha=.09)
ax.set(xlabel='Speed (km/h)',ylabel='Fuel consumption (L/100 km)',xlim=(0,160),ylim=(0,14))
ax.legend(frameon=False,fontsize=8,loc='upper center')
fig.tight_layout();fig.savefig(ROOT/'figures/fuel_curve.pdf');fig.savefig(ROOT/'figures/fuel_curve.png',dpi=180);plt.close(fig)

p=rows['coast'];vd=np.linspace(vo,vi,400);va=np.linspace(vi,vo,400)
xd=p['x12']+cumulative_trapezoid(-vd/R(vd),vd,initial=0)
xa=L+cumulative_trapezoid(va/F(va),va,initial=0)
segments=[(np.linspace(X0,p['x12'],100),np.full(100,vo),np.full(100,ucr(vo))), (xd,vd,np.zeros(len(vd))), (np.linspace(0,L,150),np.full(150,vi),np.full(150,ucr(vi))), (xa,va,np.ones(len(va))), (np.linspace(p['x45'],X1,100),np.full(100,vo),np.full(100,ucr(vo)))]
fig,(ax1,ax2)=plt.subplots(2,1,figsize=(7.1,4.9),sharex=True,gridspec_kw={'height_ratios':[1.1,1]})
for ax in [ax1,ax2]:
 ax.axvspan(0,L,color='#dce8e8',alpha=.8)
 for x in [p['x12'],0,L,p['x45']]:ax.axvline(x,color='#687580',ls='--',lw=.7)
for i,(x,v,u) in enumerate(segments):
 ax1.plot(x,v*3.6,color='#126782',lw=2)
 reg='in' if i==2 else 'out'
 # Endpoints of transitions have measure zero; use the arc's one-sided cost.
 density=(ct+cf*(qi+Gamma(v)*u)+phi(v,reg))/v*1000
 ax2.plot(x,density,color='#a44938',lw=1.7)
ax1.plot([X0,0,0,L,L,X1],[130,130,80,80,130,130],color='#777777',ls=':',lw=1,label='Posted limit')
ax1.set(ylabel='Speed (km/h)',ylim=(75,136));ax1.legend(frameon=False,loc='lower left')
ax1.text(L/2,132,'Reduced-limit zone',ha='center',fontsize=9)
ax2.set(ylabel='Instantaneous cost (CZK/km)',xlabel='Position (m)',xlim=(X0,X1),ylim=(0,13))
fig.tight_layout();fig.savefig(ROOT/'figures/speed_profile.pdf');fig.savefig(ROOT/'figures/speed_profile.png',dpi=180);plt.close(fig)
