import numpy as np
from scipy.integrate import solve_ivp,quad
from numpy.polynomial import Polynomial as Poly

def make(alpha=.1,verbose=True):
 Z=Poly([0,1]); S=1-Z*Z
 F=Poly([1+62*alpha/27,0,-8*alpha/9,0,-10*alpha/9,0,-35*alpha/27,alpha])
 D=2*F-.5*S*Poly(F.deriv().coef[1:]); A=(35*alpha/2)*Z**3*(1+2*Z-2*Z*Z)
 N=Z*A
 roots=D.roots(); residues=N(roots)/D.deriv()(roots)
 def j(z):return np.exp(np.sum(residues*np.log(1-np.asarray(z)[...,None]/roots),axis=-1)).real
 def aa(s):
  z=np.sqrt(np.maximum(1-s,0))
  a=-A(z)/(2*D(z))
  az=-(A.deriv()(z)*D(z)-A(z)*D.deriv()(z))/(2*D(z)**2)
  ass=np.divide(-az,2*z,out=np.zeros_like(z),where=z!=0)
  return a,ass
 a0,ap0=aa(0)
 def fun(s,y):
  a,ap=aa(s); c=a*a-ap
  u,up,h,hp=y
  return [up,(-(5-6*s*a)*up+(6*a-4*s*c)*u-4*c)/(2*s), hp, (-(5-6*s*a)*hp+(6*a-4*s*c)*h)/(2*s)]
 eps=1e-7; pu=-4*(a0*a0-ap0)/5; ph=6*a0/5
 sol=solve_ivp(fun,(eps,1),[pu*eps,pu,1+ph*eps,ph],rtol=2e-12,atol=2e-13,method='DOP853',dense_output=True)
 u,up,h,hp=sol.y[:,-1]; c=-(3*u+2*up)/(3*h+2*hp)
 def fields(s):
  u,up,h,hp=sol.sol(s); u=u+c*h; up=up+c*hp
  a,ap=aa(s); cc=a*a-ap
  upp=(-(5-6*s*a)*up+(6*a-4*s*cc)*u-4*cc)/(2*s)
  f=1+s*u; fp=u+s*up
  E=-3*u-2*s*up; P=-E-4*f*a
  Es=-5*up-2*s*upp; Ps=-(E+P)*(fp/(2*f)-a)
  cs=np.divide(Ps,Es,out=np.zeros_like(Ps),where=(Es!=0)&(np.asarray(s)<1))
  return f,fp,E,P,cs,Es,Ps
 fR=fields(1)[0]; FR=F(0); norm=1+26*alpha/9
 Q=quad(lambda r:F(np.sqrt(1-r*r))/j(np.sqrt(1-r*r)),0,1,epsabs=2e-13,epsrel=2e-13)[0]
 fRquad=(FR+2*Q)/(3*FR+16*alpha/9)
 ss=np.linspace(1e-5,.9999,10001); fv,fpv,E,P,cs,Es,Ps=fields(ss)
 upc=(6*a0*c-4*(a0*a0-ap0))/5
 Ec=-3*c; Pc=3*c-4*a0; csc=4*a0*(c/2-a0)/(-5*upc)
 if verbose:
  print('alpha',alpha,'C',(1-fR)/2,'fR error vs quadrature',fR-fRquad,'j0',j(1))
  print('center f2,E,P,cs',c,Ec,Pc,csc,'minf',fv.min(),'cs range',cs.min(),cs.max(),'monotonic violation',np.max(np.diff(cs)))
  print('J/OmegaR3',8*alpha/(27*norm),'min E,P',E.min(),P.min(),'max Es',Es.max())
 return F,D,A,N,j,aa,fields,fR,sol,c
if __name__=='__main__':
 for alpha in [.05,.1,.2]:make(alpha)
