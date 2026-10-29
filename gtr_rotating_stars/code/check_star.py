import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from math import log,sqrt,pi

def eos(x):
 e=x**3+1.5*x**5+x**6; p=x**5+x**6
 mu=1+2.5*x*x+2*x**3; mux=5*x+6*x*x
 cs=x*x*(5+6*x)/(3+7.5*x*x+6*x**3)
 # (e+p)/cs2 with removable zero eliminated
 ec=x*(1+2.5*x*x+2*x**3)*(3+7.5*x*x+6*x**3)/(5+6*x)
 return e,p,mu,mux,cs,ec

def build(s0):
 xc=1-19*s0/24
 def bg(x,y):
  s,q,w,d=y; e,p,mu,mux,cs,ec=eos(x)
  f=1-2*s*q; G=(q+p)/f; a=-(e+p)/(2*f); sx=-2*mux/(mu*G)
  return sx*np.array([1,(e-3*q)/(2*s),d,-((5+2*s*a)*d+4*a*w)/(2*s)])
 b=solve_ivp(bg,(xc,0),[s0,7/6-209*s0/80,1+4.4*s0,4.4],rtol=2e-11,atol=1e-12,method='DOP853',dense_output=True)
 sR,qR,wR,dR=b.y[:,-1]; R=sqrt(sR); M=R**3*qR; fR=1-2*M/R
 lam=1/(wR+2*sR*dR/3); J=R**5*lam*dR/3
 j0=5.5**2/fR; p1=j0*lam**2/3; V0=(2/3)*5.5*j0*lam**2
 z1=(8.25*p1+V0)/5; G0=19/6
 def rot(x,y):
  s,q,w,d=b.sol(x); w*=lam; d*=lam
  e,p,mu,mux,cs,ec=eos(x); f=1-2*s*q; G=(q+p)/f; a=-(e+p)/(2*f)
  sx=-2*mux/(mu*G); E=mu*mu/fR; j2=f*E
  U=(2/3)*s*j2*d*d; V=-(4/3)*j2*a*w*w
  Tp=E*((1-s*G)*w*w+2*s*w*d)/3
  z,p0,H,B,Hh,Bh=y
  zs=(ec*p0+s*(U/2+V)-3*z)/(2*s)
  ps=-(1+2*s*p)*z/(2*f*f)+a*p0+s*U/(4*f)+Tp
  C=np.array([[-2-2*s*G+(e+p-2*q)/(f*G),-2/(f*G)],[-2*G,-4]])
  D=np.array([s*s*G*(U+V)+(V-U)/(2*f*G),(1+s*G)*(U+V)])
  return sx*np.r_[zs,ps,(C@np.array([H,B])+D)/(2*s),C@np.array([Hh,Bh])/(2*s)]
 rr=solve_ivp(rot,(xc,0),[z1*s0,p1*s0,0,V0/4,1,-G0/2],rtol=2e-10,atol=1e-11,method='DOP853',dense_output=True)
 z,p0,H,B,Hh,Bh=rr.y[:,-1]
 t=R/M-1; L=log((t+1)/(t-1))
 Q22=1.5*(t*t-1)*L-(3*t**3-5*t)/(t*t-1)
 VQ=2*((3*t*t-2)/(t*t-1)-1.5*t*L)
 hp=J*J*(1/(M*R**3)+1/R**4); vp=-J*J/R**4
 mat=np.array([[sR*Hh,-Q22],[sR*sR*Bh,-VQ]])
 c,K=np.linalg.solve(mat,[hp-sR*H,vp-sR*sR*B])
 h2=sR*(H+c*Hh); v2=sR*sR*(B+c*Bh)
 dm=R**3*z+J*J/R**3; TR=sR/fR*(lam*wR)**2/3
 p2=-h2-TR; grav=M/(R*(R-2*M)); xi0=p0/grav; xi2=p2/grav
 Q=-(J*J/M+8*K*M**3/5)
 hc=-dm/(R-2*M)+J*J/(R**3*(R-2*M))+p0-TR
 out=dict(R=R,M=M,compactness=M/R,J_per_Omega=J,deltaM_per_Omega2=dm,K_per_Omega2=K,Q_per_Omega2=Q,xi0_per_Omega2=xi0,xi2_per_Omega2=xi2,shape_per_Omega2=xi2/R+v2-h2,det=np.linalg.det(mat),c=c,hc=hc)
 def local(r):
  x=brentq(lambda x:b.sol(x)[0]-r*r,0,xc,xtol=1e-14)
  s,q,w,d=b.sol(x); z,p0,H,B,Hh,Bh=rr.sol(x); e,p,mu,*_=eos(x)
  return x,q,lam*w,z,p0,s*(H+c*Hh),s*s*(B+c*Bh),mu,e,p
 def metric(r,th,om):
  x,q,w,z,p0,h2,v2,mu,e,p=local(r); s=r*r; f=1-2*s*q; E=mu*mu/fR
  # quadratic source functions
  # get W_s
  d=lam*b.sol(x)[3]; a=-(e+p)/(2*f); j2=f*E
  U=(2/3)*s*j2*d*d; V=-(4/3)*j2*a*w*w
  m2=r*f*(-h2+s*s*(U+V)); h0=hc-p0+s*E*w*w/3
  P2=(3*np.cos(th)**2-1)/2; omega=om*(1-w)
  g=np.zeros((4,4)); g[0,0]=-(1+2*om**2*(h0+h2*P2))/E+s*np.sin(th)**2*omega**2
  g[0,3]=g[3,0]=-s*np.sin(th)**2*omega
  g[1,1]=(1+2*om**2*(r**3*z+m2*P2)/(r*f))/f
  g[2,2]=s*(1+2*om**2*(v2-h2)*P2); g[3,3]=g[2,2]*np.sin(th)**2
  return g
 return out,metric,local
if __name__=='__main__':
 for s0 in [1e-6,1e-8]:
  out,metric,local=build(s0)
  print('s0',s0,out)

def einstein(metric,r,theta,om,step=1e-3):
 g=metric(r,theta,om); inv=np.linalg.inv(g)
 dg=np.zeros((4,4,4)); ddg=np.zeros((4,4,4,4))
 grid=np.arange(-2,3); c1=np.array([1,-8,0,8,-1])/12; c2=np.array([-1,16,-30,16,-1])/12
 for k in [1,2]:
  vals=np.array([metric(r+i*step*(k==1),theta+i*step*(k==2),om) for i in grid])
  dg[:,:,k]=np.einsum('i,ijk->jk',c1,vals)/step
  ddg[:,:,k,k]=np.einsum('i,ijk->jk',c2,vals)/step**2
 for i in range(5):
  for j in range(5):
   ddg[:,:,1,2]+=c1[i]*c1[j]*metric(r+grid[i]*step,theta+grid[j]*step,om)/step**2
 ddg[:,:,2,1]=ddg[:,:,1,2]
 gam=np.zeros((4,4,4)); dgam=np.zeros((4,4,4,4))
 for a in range(4):
  for b in range(4):
   for c in range(4):
    for d in range(4):
     z=dg[d,c,b]+dg[d,b,c]-dg[b,c,d]
     gam[a,b,c]+=.5*inv[a,d]*z
     for k in range(4):
      di=-(inv@dg[:,:,k]@inv)[a,d]
      dgam[a,b,c,k]+=.5*(di*z+inv[a,d]*(ddg[d,c,b,k]+ddg[d,b,c,k]-ddg[b,c,d,k]))
 Ric=np.zeros((4,4))
 for a in range(4):
  for b in range(4):
   for c in range(4):
    Ric[a,b]+=dgam[c,a,b,c]-dgam[c,a,c,b]
    for d in range(4):
     Ric[a,b]+=gam[c,a,b]*gam[d,c,d]-gam[d,a,c]*gam[c,b,d]
 return Ric-.5*g*np.sum(inv*Ric)

if __name__=='__main__':
 R=out['R']; r=.5*R; th=1.1
 x,q,w,z,p0,h2,v2,mu,e,p=local(r); fR=1-2*out['compactness']; E=mu**2/fR
 cs=eos(x)[4]; p2=-h2-r*r*E*w*w/3; basepert=p0+p2*(3*np.cos(th)**2-1)/2
 for om in [0,.04,.02,.01]:
  g=metric(r,th,om); u=np.array([1.,0.,0.,om]); u/=sqrt(-u@g@u); uc=g@u
  pt=(p+(e+p)*om**2*basepert)/(4*pi)
  et=(e+(e+p)/cs*om**2*basepert)/(4*pi)
  T=(et+pt)*np.outer(uc,uc)+pt*g
  res=einstein(metric,r,th,om)-8*pi*T
  print('Einstein residual',om,'diagonal max',np.max(np.abs(np.diag(res))),'tphi',res[0,3],'rtheta',res[1,2])
