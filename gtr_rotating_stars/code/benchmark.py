"""External rotation benchmarks.

The archived comparison integrator is used unmodified. Its Love result is
discarded for incompressible stars (a surface-density jump is required there).
"""
import sys, json, re
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

import rewrite_comparison as rc
from paths import DATA as OUT, INPUTS
original_solve_ivp=rc.solve_ivp
def checked_solve_ivp(*args,**kwargs):
    result=original_solve_ivp(*args,**kwargs)
    if not result.success:
        raise RuntimeError(result.message)
    return result
rc.solve_ivp=checked_solve_ivp

def uniform(C):
    b = np.sqrt(1-2*C)
    def fields(s):
        f = 1-2*C*s; a=np.sqrt(f); lapse=(3*b-a)/2
        E=6*C; P=E*(a-b)/(3*b-a); g=(C+P/2)/f
        return f,E,P,np.inf,0.,g,C,-(E+P)/(4*f),a/lapse
    f,E,P,v,es,g,q,aj,j=fields(0.)
    return C,fields,(E,P,v,j),1.

original_eos=rc.eos
def extended_eos(h,kind):
    if kind == 'energy_polytrope':
        E=np.expm1(np.maximum(h,0)/2)
        return E,E*E,2*E
    return original_eos(h,kind)
rc.eos=extended_eos

def energy_poly(hc):
    C,fields,central,R=rc.tov(hc,'energy_polytrope')
    def fixed_fields(s):
        vals=list(fields(s))
        if np.ndim(s)==0 and s>=1:
            vals[4]=-R*R*vals[5]/4
        return tuple(vals)
    return C,fixed_fields,central,R

def rotation(bg,s0=1e-7):
    result=rc.response(bg,s0)
    return {k:v for k,v in result.items() if k not in ('k2','cs2_c')}

def transfer_checks():
    from drag_regular import make
    from technical_family import center_data
    refs=json.loads((OUT/'family_results.json').read_text());rows=[]
    for alpha in [.1,.3]:
        F,D,A,N,j,aa,fields,fb,sol,c=make(alpha,False);C=(1-fb)/2
        Ec,Pc,v,*_=center_data(alpha,c)
        def fld(s):
            f,fs,E,P,v,Es,Ps=fields(s);g=fs/f-2*aa(s)[0];q=(1-f)/(2*s)
            return f,E,P,v,Es,g,q,aa(s)[0],j((1-s)**.5)
        out=rotation((C,fld,(Ec,Pc,v,j(1)),1.))
        ref=next(r for r in refs if r['alpha']==alpha)
        rows.append(dict(alpha=alpha,benchmark_solver=out,I_relative=out['I_R3']/ref['I']-1,Q_relative=out['Q_coeff']/ref['Q']-1))
    (OUT/'benchmark_transfer.json').write_text(json.dumps(rows,indent=2))

def calibration_checks():
    checks=[]
    for val,pos in [(.01,.3329),(.1,.7924),(.2,.9307),(.5,1.1136)]:
        pred=.15+np.log(val/.004)/np.log(150)
        checks.append(dict(axis='C',value=val,EPS_position=pos,predicted_position=pred,difference=pred-pos))
    for name,lo,hi,ticks in [('I',3,1000,[(10,.2951),(100,.5725),(1000,.8500)]),('Q',.7,300,[(1,.1912),(10,.4572),(100,.7231)])]:
        for val,pos in ticks:
            pred=.15+.7*np.log(val/lo)/np.log(hi/lo)
            checks.append(dict(axis=name,value=val,EPS_position=pos,predicted_position=pred,difference=pred-pos))
    assert max(abs(x['difference']) for x in checks)<5.01e-5
    (OUT/'EPS_axis_calibration.json').write_text(json.dumps(checks,indent=2))

def eps_curve(filename):
    if not filename.exists():
        vertices=json.loads((OUT/'published_vertices.json').read_text())
        key='I' if filename.name.startswith('I') else 'Q'
        xy=np.array([[v['source_x'],v['source_y_'+key]] for v in vertices])
    else:
        text=filename.read_text()
        # Legend text precedes its line sample: Color10 is Shen; Color11 is n=1.
        part=re.split(r'Color\d+ SC',text.split('Color11 SC',1)[1],maxsplit=1)[0]
        xy=np.array([(float(x),float(y)) for x,y in
            re.findall(r'^([\d.]+) ([\d.]+) [ml]$',part,re.M)])
        assert len(xy)>20
    # Grace world window recovered from all visible tick positions:
    # x: .004 .. .6, viewport .15 .. 1.15 (log scale).
    # y: 3 .. 1000 for I; .7 .. 300 for Q, viewport .15 .. .85.
    lo,hi=(3.,1000.) if filename.name.startswith('I') else (.7,300.)
    C=.004*np.exp((xy[:,0]-.15)*np.log(150.))
    val=lo*np.exp((xy[:,1]-.15)/.7*np.log(hi/lo))
    return xy,C,val,lo,hi

def main():
    targets=json.loads((INPUTS/'benchmarks_stage1.json').read_text())
    rows=[]
    for t in targets['rows']:
        r=rotation(uniform(t['C']))
        rr=rotation(uniform(t['C']),1e-8)
        ibar=r['I_R3']/r['C']**3
        row=dict(**t,calculated_Ibar=ibar,calculated_Qbar=r['Qbar'],
            I_relative_error=ibar/t['Ibar']-1,Q_relative_error=r['Qbar']/t['Qbar']-1,
            refinement_I_relative=rr['I_R3']/r['I_R3']-1,
            refinement_Q_relative=rr['Qbar']/r['Qbar']-1)
        print('CM',json.dumps(row),flush=True);rows.append(row)
    (OUT/'benchmark_uniform.json').write_text(json.dumps(rows,indent=2))
    # Use common source vertices, avoiding interpolation of the published curve.
    src=INPUTS/'yagi_source'
    xi,ci,vi,ilo,ihi=eps_curve(src/'IbarC-PRD.eps')
    xq,cq,vq,qlo,qhi=eps_curve(src/'QbarC-PRD.eps')
    selected=[]
    for target in [.08,.15,.2]:
        k=int(np.argmin(abs(cq-target)));C=float(cq[k]);xp=float(xq[k,0])
        ki=int(np.argmin(abs(xi[:,0]-xp)));assert xi[ki,0]==xp
        hc=brentq(lambda h:energy_poly(h)[0]-C,.001,.8,xtol=1e-12)
        bg=energy_poly(hc);r=rotation(bg);rr=rotation(bg,1e-8)
        # EPS rounding ±.00005, with conservative additional calibration
        # uncertainty .00005 per coordinate: use ±.0001 for comparison bounds.
        dx=1e-4;Clo=C*np.exp(-dx*np.log(150));Chi=C*np.exp(dx*np.log(150))
        def solve_at(c):
            h=brentq(lambda h:energy_poly(h)[0]-c,.001,.8,xtol=1e-12)
            return rotation(energy_poly(h))
        rl,rh=solve_at(Clo),solve_at(Chi)
        row=dict(C=C,hc=hc,central_sound_speed_squared=bg[2][2],
            source_x=xp,source_y_I=float(xi[ki,1]),source_y_Q=float(xq[k,1]),
            published_Ibar=float(vi[ki]),published_Qbar=float(vq[k]),
            calculated_Ibar=r['I_R3']/C**3,calculated_Qbar=r['Qbar'],
            refinement_I_relative=rr['I_R3']/r['I_R3']-1,
            refinement_Q_relative=rr['Qbar']/r['Qbar']-1)
        for key,lo,hi in [('I',ilo,ihi),('Q',qlo,qhi)]:
            pub=row['published_'+key+'bar'];calc=row['calculated_'+key+'bar']
            a,b=(rl['I_R3']/Clo**3,rh['I_R3']/Chi**3) if key=='I' else (rl['Qbar'],rh['Qbar'])
            uncertainty=max(abs(a/calc-1),abs(b/calc-1))+np.expm1(dx/.7*np.log(hi/lo))
            row[key+'_relative_error']=calc/pub-1
            row[key+'_digitization_relative_allowance']=float(uncertainty)
            row[key+'_within_digitization_allowance']=bool(abs(calc/pub-1)<uncertainty)
        print('YY',json.dumps(row),flush=True);selected.append(row)
    (OUT/'benchmark_yagi_yunes.json').write_text(json.dumps(selected,indent=2))
    (OUT/'published_vertices.json').write_text(json.dumps([
        {k:v[k] for k in ('source_x','source_y_I','source_y_Q')} for v in selected],indent=2))
    provenance_path=OUT/'benchmark_provenance.json'
    previous=json.loads(provenance_path.read_text()) if provenance_path.exists() else {}
    provenance_path.write_text(json.dumps({
        'source_hashes':previous.get('source_hashes',{}),
        'figure_numbers':'Yagi-Yunes PRD paper: Figs. 6 (I) and 7 (Q), right-hand panels',
        'uniform_source':'Chandrasekhar & Miller 1974, Table I, p.73; doi:10.1093/mnras/167.1.63',
        'polytrope_source':'Yagi & Yunes 2013, arXiv:1303.1528, original source IbarC-PRD.eps and QbarC-PRD.eps, Color11 n=1 curves (legend text precedes its line sample)',
        'polytrope_EOS':'P=K epsilon^2; scaled Pi=E^2, h=2 log(1+E)',
        'EPS_coordinate_allowance':1e-4,
        'caveat':'The digitization allowance is a conservative graphical resolution estimate, not a rigorous interval certificate or an error bound for the published ODE calculation.'},indent=2))

if __name__=='__main__':
    calibration_checks();main();transfer_checks()
