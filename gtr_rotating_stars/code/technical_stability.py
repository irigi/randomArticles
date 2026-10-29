"""Fixed-EOS sequences, without extrapolation above the reference central density.

Turning-point evidence is not a direct radial-mode eigenvalue calculation.
The EOS is kept fixed for every sequence, unlike the generating alpha-family.
"""
import json,csv
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from drag_regular import make
from technical_family import center_data,OUT

def eos_table(alpha,npoints=6001):
 bg=make(alpha,verbose=False)
 F,D,A,N,j,aa,fields,fR,sol,c=bg
 Ec,Pc,cs,css,Es=center_data(alpha,c); Gc=c-2*aa(0)[0]
 def source(z,y):
  if z>=1:return [(Ec+Pc)*Gc]
  s=1-z*z; f,fp,*_=fields(s); a=aa(s)[0]
  return [z*(-4*f*a)*(fp/f-2*a)]
 integ=solve_ivp(source,(0,1),[0.],rtol=2e-12,atol=1e-15,method='DOP853',dense_output=True)
 z=np.linspace(0,1,npoints); s=1-z*z
 f,fp,E,P,*_=fields(s[1:-1]); a=aa(s[1:-1])[0]
 P=integ.sol(z[1:-1])[0]; E=-4*f*a-P
 ent=np.r_[0,.5*np.log(fR/f)+np.log(j(z[1:-1])),.5*np.log(fR)+np.log(j(1))]
 E=np.r_[0,E,Ec]; P=np.r_[0,P,Pc]
 # Linear leading powers at the surface remove the fractional-power endpoint.
 et=PchipInterpolator(ent,E**(2/3),extrapolate=False)
 pt=PchipInterpolator(ent,P**(2/5),extrapolate=False)
 def eos(h):
  h=np.clip(h,0,ent[-1]); return et(h)**1.5,pt(h)**2.5
 return eos,float(ent[-1]),bg,float(integ.y[0,-1]-Pc)

def tov(eos,hc):
 E,P=eos(hc); q=E/6
 dh=hc*1e-7; hstart=hc-dh
 s=2*dh/(q+P/2); r=np.sqrt(s); m=q*r**3
 def rhs(h,y):
  r,m=y; E,P=eos(h)
  dr=-r*(r-2*m)/(m+.5*r**3*P)
  return [dr,.5*r*r*E*dr]
 sol=solve_ivp(rhs,(hstart,0),[r,m],rtol=2e-10,atol=2e-12,method='DOP853')
 assert sol.success
 R,M=sol.y[:,-1]
 return float(R),float(M)

def sequence(alpha,npoints=6001):
 eos,hc,bg,perror=eos_table(alpha,npoints)
 frac=np.r_[np.geomspace(.002,.1,12),np.linspace(.125,.95,34),.975,.99,.995,1.]
 rows=[]
 for f in frac:
  E,P=eos(f*hc); R,M=tov(eos,f*hc)
  rows.append(dict(alpha=alpha,hc_fraction=float(f),E_center=float(E),P_center=float(P),R=R,M=M,C=M/R))
 masses=np.array([x['M'] for x in rows]); dmdh=np.diff(masses)/np.diff(frac*hc)
 # One-sided second-order derivative entirely inside the reconstructed EOS.
 slopes={}
 for step in [.01,.005,.0025]:
  ms=[tov(eos,hc*(1-k*step))[1] for k in [0,1,2]]
  slope=(3*ms[0]-4*ms[1]+ms[2])/(2*step*hc)
  Ec,Pc,v,*_=center_data(alpha,bg[-1])
  slopes[str(step)]=slope*v/(Ec+Pc)  # dM/d(8 pi epsilon_c)
 imax=int(np.argmax(masses))
 summary=dict(alpha=alpha,hc=hc,reference_R_error=rows[-1]['R']-1,
  reference_M_error=rows[-1]['M']-(1-bg[-3])/2,pressure_integral_error=perror,
  dM_dE_reference=slopes,min_sampled_dM_dh=float(dmdh.min()),
  mass_increasing_on_sample=bool(np.all(dmdh>0)),
  maximum_sample=rows[imax],npoints=npoints)
 return rows,summary

if __name__=='__main__':
 rows=[]; summaries=[]
 for a in [.1,.2,.3,.315083901692963,1.,8.]:
  rr,ss=sequence(a);rows+=rr;summaries.append(ss);print(json.dumps(ss),flush=True)
 with (OUT/'fixed_eos_sequences.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 (OUT/'fixed_eos_stability.json').write_text(json.dumps(summaries,indent=2))
