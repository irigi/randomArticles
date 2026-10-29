"""Finite-difference Einstein residuals in the intended parameter range."""
import json
import numpy as np
from drag_second import build
from check_star import einstein
from technical_family import OUT

rows=[]
for alpha in [.1,.3]:
 out,metric,local=build(1e-8,alpha=alpha,verbose=False)
 r=.55;theta=1.1
 f,E,P,cs,w,li,z,p,h2,v2,m2,h0=local(r)
 dp=p+(-h2-r*r*li*w*w/3)*(3*np.cos(theta)**2-1)/2
 for om in [0,.04,.02,.01]:
  g=metric(r,theta,om);u=np.array([1.,0.,0.,om]);u/=np.sqrt(-u@g@u);uc=g@u
  pt=P+(E+P)*om**2*dp;et=E+(E+P)/cs*om**2*dp
  residual=einstein(metric,r,theta,om,step=4e-4)-(et+pt)*np.outer(uc,uc)-pt*g
  row=dict(alpha=alpha,r=r,theta=theta,OmegaR=om,
    diagonal_residual=float(np.max(abs(np.diag(residual)))),
    tphi_residual=float(residual[0,3]),rtheta_residual=float(residual[1,2]))
  rows.append(row);print(json.dumps(row))
(OUT/'einstein_residuals.json').write_text(json.dumps(rows,indent=2))
