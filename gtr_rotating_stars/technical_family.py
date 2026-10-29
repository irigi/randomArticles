"""Technical review calculations. No article prose is generated here.

Dimensionless units R=1, G=c=1; E=8*pi*epsilon, P=8*pi*pressure.
Results are floating-point evidence, not validated interval enclosures.
"""
import json,csv
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.optimize import minimize_scalar,brentq
from drag_regular import make
from drag_second import build

OUT=Path('technical_review');OUT.mkdir(exist_ok=True)

def center_data(alpha,c):
 a0=-35*alpha/8; a1=35*alpha/16+735*alpha**2/32
 a2=385*alpha/64-1085*alpha**2/64-15435*alpha**3/128
 u1=(6*a0*c-4*(a0*a0-a1))/5
 u2=(12*a0*u1+(10*a1-4*a0*a0)*c-8*a0*a1+8*a2)/14
 E=-3*c; P=3*c-4*a0; E1=-5*u1; E2=-7*u2
 P1=5*u1-4*(a1+c*a0); P2=7*u2-4*(a2+c*a1+u1*a0)
 cs=P1/E1; css=2*(P2*E1-E2*P1)/E1**2
 return E,P,cs,css,E1

def love(bg):
 F,D,A,N,j,aa,fields,fR,sol,c=bg
 alpha=float(A.coef[3])*2/35
 Ec,Pc,v,_,Es=center_data(alpha,c); Gc=c-2*aa(0)[0]
 chi=-2*Es/Gc
 y1=-(Ec+33*Pc+3*chi)/42
 s0=1e-7
 def rhs(s,y):
  f,fp,E,P,cs,Es,Ps=fields(s); G=fp/f-2*aa(s)[0]
  Ftid=(1-s*(E-P)/2)/f
  Qtid=(s*(5*E+9*P-2*Es/G)/2-6)/f-4*s*s*G*G
  return -(y*y+y*Ftid+Qtid)/(2*s)
 soly=solve_ivp(rhs,(s0,1),[2+y1*s0],rtol=1e-11,atol=1e-12,method='DOP853')
 assert soly.success
 y=float(soly.y[0,-1]); C=(1-fR)/2
 num=8*C**5*(1-2*C)**2*(2+2*C*(y-1)-y)/5
 den=(2*C*(6-3*y+3*C*(5*y-8))+
      4*C**3*(13-11*y+C*(3*y-2)+2*C*C*(1+y))+
      3*(1-2*C)**2*(2-y+2*C*(y-1))*np.log1p(-2*C))
 k2=num/den
 return dict(yR=y,k2=k2,Lambda=2*k2/(3*C**5))

def summary(alpha,rotation=True):
 if rotation:
  rot,metric,local,details=build(1e-8,alpha=alpha,verbose=False,details=True)
  bg=details['background']
 else:bg=make(alpha,verbose=False)
 F,D,A,N,j,aa,fields,fR,sol,c=bg
 Ec,Pc,v,css,Es0=center_data(alpha,c)
 grid=np.linspace(1e-5,.99999,3001); vals=fields(grid); cs=vals[4]
 imax=np.argmax(cs)
 if imax==0 and css<=0:smax=0; vmax=v
 else:
  lo=grid[max(0,imax-2)]; hi=grid[min(len(grid)-1,imax+2)]
  opt=minimize_scalar(lambda s:-fields(s)[4],bounds=(lo,hi),method='bounded',options={'xatol':1e-12})
  smax=float(opt.x); vmax=float(-opt.fun)
 # Near the surface P is a subtraction of O(z^3) terms leaving O(z^5).
 # Restore it from the hydrostatic integral instead of dividing cancellation noise.
 pressure=vals[3].copy()
 for k in np.flatnonzero(grid>.98):
  z=np.sqrt(1-grid[k])
  def pressure_z(t):
   s=1-t*t; f,fp,*_=fields(s); a=aa(s)[0]
   return t*(-4*f*a)*(fp/f-2*a)
  pressure[k]=quad(pressure_z,0,z,epsabs=1e-18,epsrel=1e-10)[0]
 gamma=(vals[2]+vals[3])*vals[4]/pressure
 result=dict(alpha=alpha,C=(1-fR)/2,cs2_center=v,cs2_max=vmax,
             cs2_center_slope=css,r_cs2_max=np.sqrt(smax),
             gamma_center=(Ec+Pc)*v/Pc,gamma_min=min(gamma.min(),5/3,(Ec+Pc)*v/Pc),
             max_P_over_E=np.max(vals[3]/vals[2]),min_f=vals[0].min(),
             monotone_sampled=bool(np.max(np.diff(cs))<=1e-9 and css<=1e-9))
 if rotation:
  result.update({key:float(rot[key]) for key in ['I','Q','deltaM','cond_raw','cond_scaled']})
  result.update(love(bg))
  _,_,_,_,Hh,Bh=details['rot'].sol(1)
  C=result['C']; g=C/fR
  result['yR_from_static_h2']=(-2*g-2*C/(fR*g))-2*Bh/(fR*g*Hh)
  result['love_independent_residual']=result['yR']-result['yR_from_static_h2']
  result['Ibar']=rot['I']/rot['M']**3
  result['Qbar']=-rot['Q']*rot['M']/rot['I']**2
  lnQ=np.log(result['Qbar']); lnL=np.log(result['Lambda'])
  Ifit=np.exp(np.polynomial.polynomial.polyval(lnQ,[1.35,.697,-.143,.0994,-.0124]))
  ILfit=np.exp(np.polynomial.polynomial.polyval(lnL,[1.47,.0817,.0149,.000287,-.0000364]))
  result['IQ_fractional_residual']=result['Ibar']/Ifit-1
  result['ILove_fractional_residual']=result['Ibar']/ILfit-1
  # Baryon rest mass and its O(Omega^2) correction; m_b cancels.
  def baryon(r,pert=False):
   if r==0:return 0.
   s=r*r; f,fp,E,P,cs,Es,Ps=fields(s)
   rho0=(E+P)*np.sqrt(f)/(j(np.sqrt(1-s))*np.sqrt(fR))
   val=r*r*rho0/(2*np.sqrt(f))
   if not pert:return val
   zeta,p0,*_=details['rot'].sol(s)
   W=F(np.sqrt(1-s))/(1+26*alpha/9)
   # p0/c_s^2 is integrable at the gas surface. Use open quadrature nodes.
   factor=s*zeta/f+s*j(np.sqrt(1-s))**2/f*W*W/3+p0/cs
   return val*factor
  mb=quad(lambda r:baryon(r),0,1,epsabs=2e-10,epsrel=2e-10,limit=150)[0]
  dmb=quad(lambda r:baryon(r,True),0,1,epsabs=2e-9,epsrel=2e-9,limit=150)[0]
  result.update(Mb=mb,deltaMb=dmb,deltaM_fixed_N=rot['I']/2,
                first_law_residual=rot['deltaM']-np.sqrt(fR)*dmb-rot['I']/2)
 return result

def main():
 rows=[]
 for alpha in [.03,.1,.2,.3,.35,1.,3.,8.,10.]:
  row=summary(alpha);rows.append(row)
  print(json.dumps(row),flush=True)
 with (OUT/'family_results.csv').open('w') as f:
  wr=csv.DictWriter(f,fieldnames=rows[0]);wr.writeheader();wr.writerows(rows)
 (OUT/'family_results.json').write_text(json.dumps(rows,indent=2))
 mono=brentq(lambda a:summary(a,False)['cs2_center_slope'],.3,.35,xtol=1e-10)
 causal=brentq(lambda a:summary(a,False)['cs2_max']-1,8,10,xtol=1e-9)
 thresholds={'monotonicity_candidate':summary(mono,False),'causality_candidate':summary(causal,False)}
 (OUT/'thresholds.json').write_text(json.dumps(thresholds,indent=2))
 print('Thresholds',json.dumps(thresholds),flush=True)

if __name__=='__main__':main()
