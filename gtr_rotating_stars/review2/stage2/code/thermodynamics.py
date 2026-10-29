"""Cancellation-resistant thermodynamics; no thermal EOS is assumed.

Matter variables E, Pi, rho0 are multiplied by 8*pi*R^2 (G=c=1).
theta = P/(m_b*n) is an ideal-gas-equivalent temperature, not a unique T.
"""
import sys,json,csv
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from drag_regular import make
from technical_family import center_data
OUT=Path(__file__).resolve().parents[1]/'data'

class Model:
    def __init__(self,alpha):
        self.alpha=alpha
        self.F,self.D,self.A,self.N,self.j,self.aa,self.fields,self.fb,self.sol,self.c=make(alpha,False)
        self.Ec,self.Pc,self.vc,*_=center_data(alpha,self.c)
        self.C=(1-self.fb)/2;self.gb=self.C/self.fb
        self.hc=np.log(np.sqrt(self.fb)*self.j(1))
        self.rhoc=(self.Ec+self.Pc)*np.exp(-self.hc)
        self.thetac=self.Pc/self.rhoc
        self.Kc=self.Pc/self.rhoc**(5/3)
        self.H3=35*alpha*self.fb/self.D(0)
        self.Ksurface=self.gb/(5*self.H3**(2/3))/self.Kc

    def geometry(self,z):
        s=1-z*z
        if s==0: f,fp=1.,self.c
        else: f,fp,*_=self.fields(s)
        a=-self.A(z)/(2*self.D(z))
        az=-(self.A.deriv()(z)*self.D(z)-self.A(z)*self.D.deriv()(z))/(2*self.D(z)**2)
        g=fp/f-2*a;H=-4*f*a
        Hz=-4*((-2*z*fp)*a+f*az)
        Pz=z*H*g;Ez=Hz-Pz
        return float(f),float(g),float(H),float(Pz),float(Ez)

    def hydro(self,z):
        # rescale t=z*u: small pressure remains well resolved by the quadrature.
        P=z*quad(lambda u:self.geometry(z*u)[3],0,1,epsabs=1e-28,epsrel=3e-12)[0]
        h=z*z*quad(lambda u:u*self.geometry(z*u)[1],0,1,epsabs=1e-14,epsrel=3e-12)[0]
        return P,h

    def at_z(self,z):
        r=np.sqrt(max(0,1-z*z))
        if z==0:
            return dict(r=1.,z=0.,E=0.,P=0.,rho0=0.,cs2=0.,Gamma=5/3,
                        theta=0.,enthalpy=1.,internal=0.,K=self.Ksurface)
        f,g,H,Pz,Ez=self.geometry(z)
        if z<.2:
            P,h=self.hydro(z)
        else:
            P=self.Pc if z==1 else float(self.fields(1-z*z)[3])
            h=float(np.log(np.sqrt(self.fb)*self.j(z)/np.sqrt(f)))
        rho=H*np.exp(-h);theta=P/rho
        v=self.vc if z==1 else Pz/Ez
        return dict(r=float(r),z=float(z),E=float(H-P),P=float(P),rho0=float(rho),
                    cs2=float(v),Gamma=float(H*v/P),theta=float(theta),
                    enthalpy=float(np.exp(h)),internal=float(np.expm1(h)-theta),
                    K=float(P/rho**(5/3)/self.Kc))

def main():
    summaries=[];surface=[];profiles=[]
    supplied=json.loads((ROOT/'review2/attachments/thermodynamics/computed_values.json').read_text())
    for alpha in [.1,.3]:
        m=Model(alpha)
        r=np.unique(np.r_[np.linspace(0,.98,501),1-np.geomspace(.02,5e-9,200),1.])
        vals=[m.at_z(float(np.sqrt(1-x*x))) for x in r]
        profiles.extend([dict(alpha=alpha,**v) for v in vals])
        opt=minimize_scalar(lambda z:m.at_z(z)['Gamma'],bounds=(.05,.2),method='bounded',options={'xatol':1e-11})
        dip=m.at_z(opt.x);center=m.at_z(1.)
        validation=[]
        for z in [1.,.9,.5,.1,.01,.001,.0001]:
            v=m.at_z(z);Pquad,hquad=m.hydro(z)
            f=m.geometry(z)[0]
            hmetric=float(np.log(np.sqrt(m.fb)*m.j(z)/np.sqrt(f)))
            validation.append(dict(z=z,pressure_relative=(v['P']/Pquad-1) if z>=.2 else None,
                                   pressure_check='independent algebraic vs hydrostatic' if z>=.2 else 'hydrostatic integral used directly; no independent pressure comparison',
                                   enthalpy_log_absolute=hmetric-hquad))
        for z in [.01,.003,.001,.0003,.0001]:
            v=m.at_z(z)
            surface.append(dict(alpha=alpha,z=z,Gamma=v['Gamma'],
                Gamma_linear_coefficient=(v['Gamma']-5/3)/z,
                theta_leading_ratio=v['theta']/(m.gb*z*z/5),
                theta_linear_coefficient=(v['theta']/(m.gb*z*z/5)-1)/z,
                internal_over_theta=v['internal']/v['theta']))
        # Thermodynamic identity checked with independent five-point differences.
        identity=[]
        for r0 in [.1,.3,.5,.7,.9,.97]:
            dr=1e-4
            vs=[m.at_z(np.sqrt(1-(r0+k*dr)**2)) for k in [-2,-1,1,2]]
            weights=np.array([1,-8,8,-1])/(12*dr)
            logrho_r=weights@np.log([v['rho0'] for v in vs])
            Er=weights@np.array([v['E'] for v in vs]);v0=m.at_z(np.sqrt(1-r0*r0))
            identity.append(abs(logrho_r-Er/(v0['E']+v0['P']))/max(1,abs(logrho_r)))
        row=dict(alpha=alpha,C=m.C,center=center,
            outer_local_Gamma_min=dip['Gamma'],r_at_outer_local_min=dip['r'],
            K_surface_over_K_center=m.Ksurface,
            supplied_theta_c_relative=m.thetac/supplied[str(alpha)]['central_theta_kBT_over_mb_c2']-1,
            supplied_outer_dip_absolute=dip['Gamma']-supplied[str(alpha)]['gamma_min'],
            max_thermodynamic_identity_scaled_residual=max(identity),hydrostatic_checks=validation)
        summaries.append(row);print(json.dumps(row),flush=True)
    (OUT/'thermodynamic_checks.json').write_text(json.dumps(summaries,indent=2))
    for filename,rows in [('thermodynamic_profiles.csv',profiles),('surface_asymptotics.csv',surface)]:
        with (OUT/filename).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

if __name__=='__main__':main()
