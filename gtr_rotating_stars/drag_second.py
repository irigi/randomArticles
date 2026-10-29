import numpy as np
from scipy.integrate import solve_ivp
from drag_regular import make
from check_star import einstein

def vacuum_basis(t):
 # Convergent 1/t series avoid cancellation in the logarithmic expressions.
 k=np.arange(1,81,dtype=float); y=1/t
 Q22=np.sum(8*k*(k+2)/((2*k+1)*(2*k+3))*y**(2*k+1))
 VQ=np.sum(4*(k-1)/(2*k+1)*y**(2*k))
 return Q22,VQ

def build(s0=1e-7,alpha=.1,verbose=True,details=False):
 F,D,A,N,j,aa,fields,fR,bg,c0=make(alpha,verbose=verbose)
 norm=1+26*alpha/9
 Fs=lambda z:alpha*(8/9+20*z*z/9+35*z**4/9-7*z**5/2)
 a0=aa(0)[0]; G0=c0-2*a0; j20=j(1)**2
 w0=1/norm; V0=-4*j20*a0*w0*w0/3
 p1=j20*w0*w0/3
 ec0=-5*(6*a0*c0-4*(a0*a0-aa(0)[1]))/5
 K0=-ec0/G0
 z1=(K0*p1+V0)/5
 def coeff(s):
  f,fp,E,P,cs,Es,Ps=fields(s); a=aa(s)[0]
  z=np.sqrt(max(1-s,0)); w=F(z)/norm; ws=Fs(z)/norm
  j2=j(z)**2; lapseinv=j2/f; G=fp/f-2*a; q=(1-f)/(2*s)
  U=2*s*j2*ws*ws/3; V=-4*j2*a*w*w/3
  Ts=lapseinv*((1-s*G)*w*w+2*s*w*ws)/3
  C=np.array([[-2-2*s*G+((E+P)/2-2*q)/(f*G),-2/(f*G)],[-2*G,-4]])
  d=np.array([s*s*G*(U+V)+(V-U)/(2*f*G),(1+s*G)*(U+V)])
  return f,E,P,Es,G,a,w,lapseinv,U,V,Ts,C,d
 def rhs(s,y):
  f,E,P,Es,G,a,w,li,U,V,Ts,C,d=coeff(s)
  z,p,H,B,Hh,Bh=y
  return np.r_[(-Es/G*p+s*(U/2+V)-3*z)/(2*s),
   -(1+s*P)*z/(2*f*f)+a*p+s*U/(4*f)+Ts,
   (C@np.array([H,B])+d)/(2*s),C@np.array([Hh,Bh])/(2*s)]
 rr=solve_ivp(rhs,(s0,1),[z1*s0,p1*s0,0,V0/4,1,-G0/2],method='DOP853',rtol=2e-11,atol=1e-12,dense_output=True)
 assert rr.success,rr.message
 M=(1-fR)/2; J=8*alpha/(27*norm)
 t=1/M-1; L=np.log((t+1)/(t-1))
 Q22,VQ=vacuum_basis(t)
 z,p,H,B,Hh,Bh=rr.sol(1)
 hp=J*J*(1/M+1); vp=-J*J
 mat=np.array([[Hh,-Q22],[Bh,-VQ]])
 scales=np.linalg.norm(mat,axis=0)
 c,K=np.linalg.solve(mat/scales,[hp-H,vp-B])/scales
 h2=H+c*Hh; v2=B+c*Bh
 dm=z+J*J; T=F(0)**2/(3*fR*norm**2)
 hc=-z/fR+p-T; xi0=p*fR/M; xi2=(-h2-T)*fR/M
 out=dict(M=M,I=J,deltaM=dm,K=K,Q=-J*J/M-8*K*M**3/5,
  xi0=xi0,xi2=xi2,k2surface=v2-h2,shape=xi2+v2-h2,c=c,hc=hc,
  surface_p0=p,surface_h2=h2,surface_v2=v2,det=np.linalg.det(mat),
  cond_raw=np.linalg.cond(mat),cond_scaled=np.linalg.cond(mat/scales))
 def local(r):
  s=r*r; f,E,P,Es,G,a,w,li,U,V,Ts,C,d=coeff(s)
  z,p,H,B,Hh,Bh=rr.sol(s)
  h2=s*(H+c*Hh); v2=s*s*(B+c*Bh)
  m2=r*f*(-h2+s*s*(U+V)); h0=hc-p+s*li*w*w/3
  return f,E,P,fields(s)[4],w,li,z,p,h2,v2,m2,h0
 def metric(r,th,om):
  f,E,P,cs,w,li,z,p,h2,v2,m2,h0=local(r); s=r*r
  P2=(3*np.cos(th)**2-1)/2; omega=om*(1-w)
  g=np.zeros((4,4)); g[0,0]=-(1+2*om**2*(h0+h2*P2))/li+s*np.sin(th)**2*omega**2
  g[0,3]=g[3,0]=-s*np.sin(th)**2*omega
  g[1,1]=(1+2*om**2*(r**3*z+m2*P2)/(r*f))/f
  g[2,2]=s*(1+2*om**2*(v2-h2)*P2); g[3,3]=g[2,2]*np.sin(th)**2
  return g
 if details:return out,metric,local,dict(rot=rr,coeff=coeff,background=(F,D,A,N,j,aa,fields,fR,bg,c0))
 return out,metric,local

if __name__=='__main__':
 for s0 in [1e-7,1e-9]:
  out,metric,local=build(s0)
  print('second order',s0,out)
 r=.55; th=1.1
 f,E,P,cs,w,li,z,p,h2,v2,m2,h0=local(r)
 dp=p+(-h2-r*r*li*w*w/3)*(3*np.cos(th)**2-1)/2
 for om in [0,.02,.01,.005]:
  g=metric(r,th,om); u=np.array([1.,0.,0.,om]); u/=np.sqrt(-u@g@u); uc=g@u
  pt=P+(E+P)*om**2*dp; et=E+(E+P)/cs*om**2*dp
  res=einstein(metric,r,th,om,step=4e-4)-(et+pt)*np.outer(uc,uc)-pt*g
  print('residual',om,'diagonal',max(abs(np.diag(res))),'tphi',res[0,3],'rtheta',res[1,2])
