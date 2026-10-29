"""Matched-compactness comparison models, in full GR.
All background matter variables are E=8*pi*epsilon and P=8*pi*pressure.
The polytrope index below is a polytropic index, not particle density.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from drag_regular import make
from drag_second import vacuum_basis

OUT=Path('rewrite');OUT.mkdir(exist_ok=True)

def eos(h,kind):
 h=np.maximum(h,0)
 if kind=='polytrope':
  rho=np.expm1(h)/2;P=rho*rho;E=rho+P;v=2*rho/(1+2*rho)
 else:
  u=-3*np.expm1(-2*h);P=u*u;E=12*u-5*u*u;v=u/(6-5*u)
 return E,P,v

def tov(hc,kind):
 Ec,Pc,vc=eos(hc,kind);g0=Ec/6+Pc/2
 r0=1e-6/np.sqrt(Ec)
 def fun(r,y):
  m,h=y;E,P,v=eos(h,kind)
  return [.5*r*r*E,-(m+.5*r**3*P)/(r*(r-2*m))]
 def event(r,y):return y[1]
 event.terminal=True;event.direction=-1
 sol=solve_ivp(fun,(r0,100/np.sqrt(Ec)),[Ec*r0**3/6,hc-g0*r0*r0/2],
               method='DOP853',rtol=2e-11,atol=1e-13,events=event,dense_output=True)
 assert len(sol.t_events[0])==1
 R=sol.t_events[0][0];M=sol.y_events[0][0,0];C=M/R;fb=1-2*C
 def fields(s):
  r=R*np.sqrt(s);m,h=sol.sol(r);E,P,v=eos(h,kind);E*=R*R;P*=R*R
  f=1-2*m/r;q=m*R*R/r**3;g=(q+P/2)/f
  ps=-(E+P)*g/2
  es=np.divide(ps,v,out=np.zeros_like(ps),where=v!=0)
  # These two EOSs have E proportional to h at the surface.
  if np.ndim(s)==0 and s>=1:
   es=-(R*R/2 if kind=='polytrope' else 72*R*R)*g/2
  a=-(E+P)/(4*f);j=np.sqrt(f/fb)*np.exp(h)
  return f,E,P,v,es,g,q,a,j
 return C,fields,(Ec*R*R,Pc*R*R,vc,1/np.sqrt(fb)*np.exp(hc)),R

def tolman(C):
 def fun(s,y):
  P,h=y;E=15*C*(1-s);q=C*(5-3*s)/2;f=1-2*s*q;g=(q+P/2)/f
  return [-(E+P)*g/2,-g/2]
 sol=solve_ivp(fun,(1,0),[0,0],rtol=2e-12,atol=2e-13,method='DOP853',dense_output=True)
 def fields(s):
  P,h=sol.sol(s);E=15*C*(1-s);q=C*(5-3*s)/2;f=1-2*s*q;g=(q+P/2)/f
  es=-15*C;v=(E+P)*g/(30*C);a=-(E+P)/(4*f);j=np.sqrt(f/(1-2*C))*np.exp(h)
  return f,E,P,v,es,g,q,a,j
 f,E,P,v,es,g,q,a,j=fields(0)
 return C,fields,(E,P,v,j),1.

def response(bg,s0=1e-7):
 C,fields,central,R=bg;Ec,Pc,vc,jc=central;g0=Ec/6+Pc/2;a0=-(Ec+Pc)/4
 def dragging(s,y):
  a=fields(s)[7];w,ws=y
  return [ws,-((5+2*s*a)*ws+4*a*w)/(2*s)]
 ws0=-4*a0/5
 drag=solve_ivp(dragging,(s0,1),[1+ws0*s0,ws0],rtol=3e-11,atol=1e-12,method='DOP853',dense_output=True)
 wR,wsR=drag.y[:,-1];norm=wR+2*wsR/3;I=wsR/(3*norm)
 def rhs(s,y):
  f,E,P,v,es,g,q,a,j=fields(s);w,ws=drag.sol(s)/norm
  U=2*s*j*j*ws*ws/3;V=-4*j*j*a*w*w/3
  mat=np.array([[-2-2*s*g+((E+P)/2-2*q)/(f*g),-2/(f*g)],[-2*g,-4]])
  d=np.array([s*s*g*(U+V)+(V-U)/(2*f*g),(1+s*g)*(U+V)])
  return np.r_[(mat@y[:2]+d)/(2*s),mat@y[2:]/(2*s)]
 V0=-4*jc*jc*a0/(3*norm*norm)
 rot=solve_ivp(rhs,(s0,1),[0,V0/4,1,-g0/2],rtol=3e-10,atol=2e-12,method='DOP853')
 H,B,Hh,Bh=rot.y[:,-1];Qv,Vv=vacuum_basis(1/C-1)
 mat=np.array([[Hh,-Qv],[Bh,-Vv]]);sc=np.linalg.norm(mat,axis=0)
 c,K=np.linalg.solve(mat/sc,[I*I*(1/C+1)-H,-I*I-B])/sc
 Q=-I*I/C-8*K*C**3/5
 def tide(s,y):
  f,E,P,v,es,g,q,a,j=fields(s)
  ft=(1-s*(E-P)/2)/f;qt=(s*(5*E+9*P-2*es/g)/2-6)/f-4*s*s*g*g
  return -(y*y+ft*y+qt)/(2*s)
 chi=(Ec+Pc)/vc;y1=-(Ec+33*Pc+3*chi)/42
 tide_sol=solve_ivp(tide,(s0,1),[2+s0*y1],rtol=3e-11,atol=1e-12,method='DOP853')
 y=tide_sol.y[0,-1]
 num=8*C**5*(1-2*C)**2*(2+2*C*(y-1)-y)/5
 den=2*C*(6-3*y+3*C*(5*y-8))+4*C**3*(13-11*y+C*(3*y-2)+2*C*C*(1+y))+3*(1-2*C)**2*(2-y+2*C*(y-1))*np.log1p(-2*C)
 k2=num/den
 return dict(C=float(C),I_R3=float(I),I_MR2=float(I/C),Q_coeff=float(Q),Qbar=float(-Q*C/I**2),k2=float(k2),cs2_c=float(vc))

if __name__=='__main__':
 rows=[]
 for alpha in [.1,.3]:
  bg=make(alpha,False);F,D,A,N,j,aa,fields,fb,sol,c=bg;C=(1-fb)/2
  for kind in ['tolman_vii','polytrope','buchdahl']:
   if kind=='tolman_vii':model=tolman(C);hc=None
   else:
    hc=brentq(lambda h:tov(h,kind)[0]-C,.002,.2 if kind=='buchdahl' else .8,xtol=2e-12)
    model=tov(hc,kind)
   row=dict(alpha_reference=alpha,model=kind,hc=hc,**response(model))
   if kind=='buchdahl':row['C_exact_residual']=row['C']-float(-np.expm1(-2*hc)/2)
   print(json.dumps(row),flush=True);rows.append(row)
 (OUT/'comparison_results.json').write_text(json.dumps(rows,indent=2))
