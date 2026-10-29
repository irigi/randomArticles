"""Numerical Rayleigh--Ritz radial pulsation spectrum of the STATIC star.

Positive Ritz eigenvalues alone are upper bounds, not rigorous stability proofs.
Convergence and independent fixed-EOS sequences supply numerical evidence.
"""
import json
import numpy as np
from numpy.polynomial.legendre import Legendre,leggauss
from scipy.linalg import eigh
from drag_regular import make
from technical_family import OUT

def spectrum(alpha,degree=18,nquad=700):
 F,D,A,N,j,aa,fields,fR,sol,c=make(alpha,verbose=False)
 x,w=leggauss(nquad);r=(1+x)/2;w=w/2;s=r*r;z=np.sqrt(1-s)
 f,fp,E,P,v,Es,Ps=fields(s); lapse=np.sqrt(f)/j(z)
 ep=-4*f*aa(s)[0]; pr=2*r*Ps
 # Common factor 8*pi cancels from the generalized eigenvalue problem.
 pref=lapse**3/(np.sqrt(f)*s)
 Pi=pref*ep*v
 Q=pref*(pr*pr/ep-4*pr/r-P*ep/f)
 W=ep*lapse/(f**1.5*s)
 basis=[];deriv=[]
 for k in range(degree):
  pol=Legendre.basis(k); t=2*s-1
  basis.append(r**3*pol(t))
  deriv.append(3*s*pol(t)+4*s*s*pol.deriv()(t))
 Z=np.array(basis);DZ=np.array(deriv)
 K=(DZ*(w*Pi))@DZ.T-(Z*(w*Q))@Z.T
 B=(Z*(w*W))@Z.T
 vals=eigh(K,B,eigvals_only=True)
 return dict(alpha=alpha,degree=degree,nquad=nquad,
             sigma2_R2=vals[:3].tolist(),mass_matrix_condition=float(np.linalg.cond(B)))

if __name__=='__main__':
 rows=[]
 for a in [.1,.2,.3,.315083901692963,1.,8.]:
  for degree,nquad in [(8,350),(16,700),(24,1000)]:
   row=spectrum(a,degree,nquad); rows.append(row);print(json.dumps(row),flush=True)
 (OUT/'radial_modes.json').write_text(json.dumps(rows,indent=2))
 from scipy.optimize import brentq
 from technical_family import summary
 from technical_stability import eos_table,tov
 alpha=brentq(lambda a:spectrum(a,16,350)['sigma2_R2'][0],.32,.5,xtol=1e-9)
 result=dict(alpha=alpha,radial=spectrum(alpha,24,700),background=summary(alpha,False))
 eos,hc,bg,pe=eos_table(alpha)
 for step in [.005,.0025,.00125]:
  masses=[tov(eos,hc*(1-k*step))[1] for k in [0,1,2]]
  result['dM_dh_'+str(step)]=(3*masses[0]-4*masses[1]+masses[2])/(2*step*hc)
 (OUT/'radial_threshold.json').write_text(json.dumps(result,indent=2))
