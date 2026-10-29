"""Independent Volterra-series evaluation and physical/matching checks.

All quantities use R=1; rotation coefficients are divided by Omega**2.
The quadrature discretization is numerical, not a rigorous interval certificate.
"""
import json
import numpy as np
from scipy.integrate import cumulative_simpson
from drag_regular import make
from drag_second import build

def volterra_check(n=4096, iterations=45):
    F,D,A,N,j,aa,fields,fR,bg,c=make(.1)
    norm=58/45
    u=np.linspace(0,1,n+1); s=np.sin(np.pi*u/2)**2
    ds=np.pi*np.sin(np.pi*u)/2
    sc=np.maximum(s,1e-9); z=np.sqrt(np.maximum(1-sc,0))
    with np.errstate(invalid='ignore'):
        f,fp,E,P,cs,Es,Ps=fields(sc)
    a=aa(sc)[0]; G=fp/f-2*a; q=(1-f)/(2*sc)
    w=F(z)/norm
    ws=(4/45+2*z*z/9+7*z**4/18-7*z**5/20)/norm
    j2=j(z)**2; U=2*sc*j2*ws*ws/3; V=-4*j2*a*w*w/3
    Ts=j2/f*((1-sc*G)*w*w+2*sc*w*ws)/3
    b=c-2*aa(0)[0]; chi0=-(-6*aa(0)[0]*c+4*(aa(0)[0]**2-aa(0)[1]))/b
    # E_s(0)=-6*a0*c+4*(a0*a0-ap0)
    V0=-4*j(1)**2*aa(0)[0]/(3*norm**2)
    Cm=np.zeros((n+1,2,2)); dm=np.zeros((n+1,2))
    Cm[:,0,0]=-3; Cm[:,0,1]=-Es/G
    Cm[:,1,0]=-sc*(1+sc*P)/f**2; Cm[:,1,1]=2*sc*a
    dm[:,0]=sc*(U/2+V); dm[:,1]=sc*sc*U/(2*f)+2*sc*Ts
    Cm0=np.array([[-3,chi0],[0,0.]])
    Cq=np.zeros_like(Cm); dq=np.zeros_like(dm)
    Cq[:,0,0]=-2-2*sc*G+((E+P)/2-2*q)/(f*G)
    Cq[:,0,1]=-2/(f*G); Cq[:,1,0]=-2*G; Cq[:,1,1]=-4
    dq[:,0]=sc*sc*G*(U+V)+(V-U)/(2*f*G)
    dq[:,1]=(1+sc*G)*(U+V)
    Cq0=np.array([[-1,-2/b],[-2*b,-4]])
    def solve(C,d,C0,nu,y0):
        L=(C-C0)/(2*sc[:,None,None])
        d0=-C0@y0
        g=(np.einsum('nij,j->ni',C-C0,y0)+d-d0)/(2*sc[:,None])
        P0=np.eye(2)+C0/nu; Pm=-C0/nu
        def integrate(v):
            i0=cumulative_simpson(v*ds[:,None],x=u,axis=0,initial=0)
            im=cumulative_simpson(v*(ds*s**(nu/2))[:,None],x=u,axis=0,initial=0)
            im[1:]/=s[1:,None]**(nu/2)
            return i0@P0.T+im@Pm.T
        term=integrate(g); total=term.copy(); history={}
        for k in range(1,iterations+1):
            term=integrate(np.einsum('nij,nj->ni',L,term)); total+=term
            if k in [5,10,20,30,45]:
                history[k]=float(np.max(np.abs(term)))
        return total+y0,history
    ym,hm=solve(Cm,dm,Cm0,3,np.zeros(2))
    yp,hp=solve(Cq,dq,Cq0,5,np.array([0,V0/4]))
    yh,hh=solve(Cq,np.zeros_like(dq),Cq0,5,np.array([1,-b/2]))
    M=(1-fR)/2; J=2/87; t=1/M-1; log=np.log((t+1)/(t-1))
    Q22=1.5*(t*t-1)*log-(3*t**3-5*t)/(t*t-1)
    VQ=2*((3*t*t-2)/(t*t-1)-1.5*t*log)
    c2,K=np.linalg.solve(np.array([[yh[-1,0],-Q22],[yh[-1,1],-VQ]]),
                        np.array([J*J*(1/M+1),-J*J])-yp[-1])
    h2,v2=yp[-1]+c2*yh[-1]
    result=dict(grid=n,deltaM=float(ym[-1,0]+J*J),K=float(K),
                Q=float(-J*J/M-8*K*M**3/5),p0R=float(ym[-1,1]),
                h2R=float(h2),v2R=float(v2),terms=dict(monopole=hm,quad_part=hp,quad_hom=hh))
    print(json.dumps(result,indent=2))
    return result

def eos_and_matching_check():
    F,D,A,N,j,aa,fields,fR,bg,c=make(.1)
    a0=-7/16; a1=287/640; a2=7973/25600
    u1=(6*a0*c-4*(a0*a0-a1))/5
    u2=(12*a0*u1+(10*a1-4*a0*a0)*c-8*a0*a1+8*a2)/14
    E1=-5*u1; E2=-7*u2
    P1=5*u1-4*(a1+c*a0); P2=7*u2-4*(a2+c*a1+u1*a0)
    v=P1/E1; slope=2*(P2*E1-E2*P1)/E1**3
    print('central EOS extension (8piR2 energy units)',dict(v=v,dv_dE=slope,
        lambda_=slope/(1-v),E_c=-3*c,P_c=3*c-4*aa(0)[0]))
    # The homogeneous scalar quadrupole operator, applied to a natural dragging ansatz.
    # v_trial=r^4*F; use analytic F derivatives and background g'.
    for r in [.25,.55,.85]:
        s=r*r; z=np.sqrt(1-s); f,fp,E,P,cs,Es,Ps=fields(s)
        a,ap=aa(s); G=fp/f-2*a; g=r*G
        u,up,h,hp=bg.sol(s); u+=c*h; up+=c*hp
        upp=(-(5-6*s*a)*up+(6*a-4*s*(a*a-ap))*u-4*(a*a-ap))/(2*s)
        fpp=2*up+s*upp; Gs=fpp/f-fp*fp/f**2-2*ap; gp=G+2*s*Gs
        q=(1-f)/(2*s); L=gp/g-2*g+((E+P)/2-2*q)/(f*g)
        Fr=-r*F.deriv()(z)/z
        Frr=-(1/z+r*r/z**3)*F.deriv()(z)+r*r/z**2*F.deriv(2)(z)
        trial=r**4*F(z); tp=4*r**3*F(z)+r**4*Fr
        tpp=12*r*r*F(z)+8*r**3*Fr+r**4*Frr
        print('trial v=r4 F scalar residual',r,tpp-L*tp-4*trial/(r*r*f))
    out,metric,local=build(1e-9)
    M=out['M']; J=out['I']; K=out['K']
    def exterior(r):
        x=r/M-1; L=np.log((x+1)/(x-1))
        Q=1.5*(x*x-1)*L-(3*x**3-5*x)/(x*x-1)
        V=2*((3*x*x-2)/(x*x-1)-1.5*x*L)
        return np.array([J*J*(1/(M*r**3)+1/r**4)+K*Q,-J*J/r**4+K*V])
    r=1.; f,E,P,cs,w,li,z,p,h,v,m2,h0=local(r)
    wr=2*(4/45)/(58/45); g=M/fR
    Ar=wr*wr/6; Br=0
    hp=(-2*g-2*M/(f*g))*h-2*v/(f*g)+(g-1/(2*f*g))*Ar
    vp=-2*g*h+(1+g)*Ar
    step=.002
    ed=(exterior(1-2*step)-8*exterior(1-step)+8*exterior(1+step)-exterior(1+2*step))/(12*step)
    print('surface derivative differences h2,v2',np.array([hp,vp])-ed)

if __name__=='__main__':
    for n in [2048,4096]:volterra_check(n)
    eos_and_matching_check()
