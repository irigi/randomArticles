"""Direct-transcription solvers used for Figs. 1-2 and value-of-information calculations.

Uniform-prior expectations use exact probability mass and exact survival-time
integrals over each time cell.  This avoids the narrow-prior quadrature bug that
occurs if midpoint density samples are used on an unaligned grid.
"""
from __future__ import annotations
import numpy as np
import casadi as ca
from pathlib import Path
from model import Model, D0, VMAX, POST, road_load


def _rk4_step(model, v, d, u, b, dt, return_stages=False):
    def deriv(vv):
        return model.rhs_v(vv, u, b)
    k1 = deriv(v)
    k2 = deriv(v + 0.5*dt*k1)
    k3 = deriv(v + 0.5*dt*k2)
    k4 = deriv(v + dt*k3)
    vn = v + dt*(k1 + 2*k2 + 2*k3 + k4)/6
    # Consistent RK4 quadrature for d_dot=-v.
    dn = d - dt*(v + 2*(v+0.5*dt*k1) + 2*(v+0.5*dt*k2) + (v+dt*k3))/6
    if return_stages:
        return vn, dn, (v + 0.5*dt*k1, v + 0.5*dt*k2, v + dt*k3)
    return vn, dn


def _known_initializations(tau, N, model, warm_starts=None):
    """Physically interpretable initializations for the fixed-switch NLP.

    All arrays are *initial guesses*, not asserted solutions.  The optimizer
    enforces exact RK4 dynamics and the red-line state constraint afterwards.
    In particular, rolling guesses place the vehicle well short of the line.
    """
    t = np.linspace(0, tau, N+1)
    dt = tau/N

    def make(label, vg, ug=None, bg=None):
        vg = np.maximum(0., np.minimum(VMAX, np.asarray(vg, dtype=float)))
        # Initial d follows initial v, and stays in the admissible domain.
        dg = D0 - np.r_[0., np.cumsum((vg[:-1]+vg[1:])*dt/2)]
        if ug is None or bg is None:
            vdot = np.diff(vg)/dt
            vr = 0.5*(vg[:-1]+vg[1:])
            desired = vdot + np.asarray(road_load(vr))
            ug = np.clip(desired/np.asarray(model.throttle_accel(vr)), 0., 1.)
            bg = np.clip(-desired/model.brake_a, 0., 1.)
        return (label, vg, np.maximum(dg, 0.), np.asarray(ug), np.asarray(bg))

    # The original single start, retained as a direct regression control.
    vg = VMAX*(1.0 - .7*t/max(tau, 1e-9))
    yield make('legacy', vg, np.zeros(N), np.full(N, .04))

    # Brake quickly to a sustainable rolling speed, hold it through most of
    # the red, and (optionally) rebuild speed approaching the green change.
    # Distinct speed levels distinguish slow rolling from stop/wait basins.
    cruise = D0/tau
    for frac in (0.65, 0.95):
        floor = min(VMAX*.8, max(.8, frac*cruise))
        ramp = max(.65, (VMAX-floor)/(model.brake_a+float(road_load(VMAX))))
        vg = floor + (VMAX-floor)*np.clip(1-t/ramp, 0., 1.)
        if frac == .95:
            # Recovery shortly before green provides a separate control basin.
            acc_time = max(.5, (VMAX-floor)/max(float(model.net_full_accel(floor)), 1.))
            vg += (VMAX-floor)*np.clip((t-(tau-acc_time))/acc_time, 0., 1.)
        yield make(f'rolling_{frac:.2f}', np.minimum(vg, VMAX))

    # Fast braking to a complete stop well before the switch; the original
    # transcription permits holding speed at zero by force-balancing throttle.
    ramp = max(.65, VMAX/(model.brake_a+float(road_load(VMAX))))
    vg = VMAX*np.clip(1.0-t/ramp, 0., 1.)
    yield make('stop_wait', vg)

    # Coast for some time, then brake; useful for late arrival without
    # assuming an immediate severe braking maneuver.
    coast_until = min(.38*tau, max(1.0, (tau-D0/VMAX)*.5))
    tail = np.clip((t-coast_until)/max(.75,ramp), 0., 1.)
    vg = VMAX*(1-tail) - np.minimum(t,coast_until)*float(road_load(VMAX))
    yield make('coast_then_stop', np.maximum(vg, 0.))

    # Nearby-switch solutions map to the same normalized time, useful both
    # for S increasing and decreasing. Their feasibility is re-enforced by NLP.
    for j, prev in enumerate(warm_starts or []):
        tv = np.asarray(prev['t'],float)
        ratio = tv/max(float(tv[-1]),1e-12)
        vg = np.interp(t/tau, ratio, np.asarray(prev['v'],float))
        ug = np.interp((t[:-1]+dt*.5)/tau,
                       (tv[:-1]+np.diff(tv)*.5)/tv[-1],np.asarray(prev['u'],float))
        bg = np.interp((t[:-1]+dt*.5)/tau,
                       (tv[:-1]+np.diff(tv)*.5)/tv[-1],np.asarray(prev['b'],float))
        yield make(f'continuation_{j}', vg, ug, bg)


def solve_known(tau, N=None, model=None, *, multistart=True,
                warm_starts=None, return_candidates=True,
                start_labels=None, max_cpu_time=12.0):
    """Best converged *local* candidate found via deterministic multistart.

    ``multistart=False`` preserves the legacy initialization for regression.
    Pass nearby solutions as ``warm_starts`` for continuation in either
    direction of switch time. Solver metadata records *all* trial outcomes.
    Lower costs are evidence of improvements, not proof of global optimality.
    """
    model = model or Model(300.0)
    if N is None:
        N = max(180, int(round(12*tau)))
    dt = tau/N
    op = ca.Opti()
    v = op.variable(N+1); d = op.variable(N+1)
    u = op.variable(N); b = op.variable(N)
    op.subject_to(v[0] == VMAX); op.subject_to(d[0] == D0)
    op.subject_to(op.bounded(0, v, VMAX)); op.subject_to(op.bounded(0, d, D0))
    op.subject_to(op.bounded(0, u, 1)); op.subject_to(op.bounded(0, b, 1))
    obj = 0
    for k in range(N):
        vn, dn, stages = _rk4_step(model, v[k], d[k], u[k], b[k], dt,
                                   return_stages=True)
        # Red-line viability also requires nonnegative intermediate RK4
        # speeds; with nonnegative speed the distance d cannot decrease
        # below its next-node value anywhere within a cell.
        for stage_v in stages:
            op.subject_to(stage_v >= 0)
        op.subject_to(v[k+1] == vn); op.subject_to(d[k+1] == dn)
        obj += dt*model.running_cost((v[k]+v[k+1])/2, u[k])
    obj += model.green_value(d[N], v[N])
    op.minimize(obj)
    op.solver('ipopt', {'expand':True,'print_time':False},
              {'print_level':0,'max_iter':3000,'tol':1e-9,'acceptable_tol':1e-8,
               'bound_relax_factor':0,'sb':'yes','mu_strategy':'adaptive',
               'warm_start_init_point':'no', 'max_cpu_time':max_cpu_time})
    t = np.linspace(0, tau, N+1)
    candidates = []
    best = None
    initializers = _known_initializations(tau,N,model,warm_starts)
    for label,vg,dg,ug,bg in initializers:
        if not multistart and label != 'legacy':
            break
        if start_labels is not None and label not in start_labels:
            continue
        op.set_initial(v,vg); op.set_initial(d,dg)
        op.set_initial(u,ug); op.set_initial(b,bg)
        try:
            sol = op.solve()
            val=lambda z:np.asarray(sol.value(z)).ravel()
            cost = float(sol.value(obj))
            V,D,U,B=val(v),val(d),val(u),val(b)
            # Independently evaluate the collocation defects and path bounds.
            # These are local feasibility checks, *not* optimality certificates.
            defect = 0.
            min_stage_speed = float(min(V))
            for k in range(N):
                vv,dd,stage = _rk4_step(model,float(V[k]),float(D[k]),
                                  float(U[k]),float(B[k]),dt,
                                  return_stages=True)
                min_stage_speed = min(min_stage_speed, *[float(z) for z in stage])
                defect=max(defect,abs(float(vv)-V[k+1]),abs(float(dd)-D[k+1]))
            feasibility = max(defect,max(0.,-min(V)) ,max(0.,-min(D)),
                              max(0.,max(V)-VMAX), max(0.,max(D)-D0),
                              max(0.,-min_stage_speed))
            if not np.isfinite(cost) or feasibility > 2e-5:
                raise ValueError(f'failed independent feasibility check: {feasibility}')
            candidate = {'t':t,'v':V,'d':D,'u':U,'b':B,'cost':cost,'T':tau}
            candidates.append({'start':label, 'success':True, 'cost':cost,
                               'min_speed':float(min(V)),
                               'min_distance':float(min(D)),
                               'min_stage_speed':min_stage_speed,
                               'max_defect':float(defect)})
            if best is None or cost < best['cost']-1e-9:
                best = candidate
                best['selected_start']=label
        except (RuntimeError, ValueError) as exc:
            candidates.append({'start':label,'success':False,'reason':str(exc)[:250]})
    if best is None:
        raise RuntimeError(f'No feasible known-time solution for S={tau}; diagnostics: {candidates}')
    if return_candidates:
        best['candidates']=candidates
    return best


def solve_known_grid(times, N=480, model=None, *, both_directions=True,
                     max_cpu_time=10., start_labels=None):
    """Warm-start sweep; re-visits points from high to low switch times.

    Returns solutions in the caller's time ordering. All times in seconds.
    The reverse sweep also runs independent multistarts, not a single warm
    start, to avoid accepting continuation artifacts as optima.
    """
    model = model or Model(300.0)
    ts = list(map(float,times))
    if len(set(ts)) != len(ts):
        raise ValueError('Switch times must be unique')
    solved = {}
    for order in (sorted(ts), sorted(ts, reverse=True)) if both_directions else (sorted(ts),):
        prev = None
        for tau in order:
            warms = [prev] if prev is not None else []
            if tau in solved:
                warms.append(solved[tau])
            result = solve_known(tau,N=N,model=model,warm_starts=warms,
                                 max_cpu_time=max_cpu_time,
                                 start_labels=(start_labels + [f'continuation_{j}' for j in range(len(warms))])
                                 if start_labels is not None else None)
            if tau not in solved or result['cost'] < solved[tau]['cost']:
                solved[tau] = result
            prev = solved[tau]
    return [solved[t] for t in ts]


def _uniform_cell_weights(start, end, lo, hi):
    """Return integral survival dt and switch probability mass for one cell."""
    before = max(0.0, min(end, lo)-start)
    a = max(start, lo); z = min(end, hi)
    inside = max(0.0, z-a)
    survival_integral = before
    if inside > 0:
        survival_integral += inside*(hi-(a+z)/2)/(hi-lo)
    mass = inside/(hi-lo)
    return survival_integral, mass


def _uniform_initializations(lo, hi, N, model, warm_starts=None):
    """Physical initial guesses for the conditional-red uniform-prior NLP.

    The guesses deliberately span qualitatively different basins: gradual
    slowing, rolling delay, stop/wait, late braking, and brake/recovery.  They
    are only starts; the NLP re-enforces dynamics and path constraints.
    Nearby prior solutions can be supplied through ``warm_starts`` and are
    rescaled onto the current horizon.
    """
    t = np.linspace(0.0, hi, N+1)
    dt = hi/N

    def make(label, vg, ug=None, bg=None):
        vg = np.maximum(0., np.minimum(VMAX, np.asarray(vg, dtype=float)))
        dg = D0 - np.r_[0., np.cumsum((vg[:-1]+vg[1:])*dt/2)]
        if ug is None or bg is None:
            vdot = np.diff(vg)/dt
            vr = 0.5*(vg[:-1]+vg[1:])
            desired = vdot + np.asarray(road_load(vr))
            ug = np.clip(desired/np.asarray(model.throttle_accel(vr)), 0., 1.)
            bg = np.clip(-desired/model.brake_a, 0., 1.)
            # Under the smooth publication convention, zero-speed waiting is
            # held by balancing road load.  In contact mode, zero controls are
            # the physically natural rest guess.
            stopped = (vg[:-1] < 1e-6) & (vg[1:] < 1e-6)
            if np.any(stopped):
                if model.rest_mode == 'contact':
                    ug[stopped] = 0.; bg[stopped] = 0.
                else:
                    ug[stopped] = np.asarray(model.cruise_throttle(0.0))
                    bg[stopped] = 0.
        return (label, vg, np.maximum(dg, 0.), np.asarray(ug), np.asarray(bg))

    # Exact legacy starting point for regression.
    vg = VMAX*(1.0-.75*t/max(hi,1e-12))
    yield make('legacy', vg, np.zeros(N), np.full(N,.04))

    # Two rolling-delay basins.  The horizon-average speed D0/hi is the scale
    # required to arrive near the line only at the latest possible switch.
    avg = D0/max(hi,1e-12)
    for frac in (0.65, 1.05):
        floor = min(.82*VMAX, max(.8, frac*avg))
        ramp = max(.65, (VMAX-floor)/(model.brake_a+float(road_load(VMAX))))
        vg = floor + (VMAX-floor)*np.clip(1-t/ramp,0.,1.)
        yield make(f'rolling_{frac:.2f}', vg)

    # Stop/wait basin, useful for broad priors where a long residual red is
    # possible.  This is intentionally distinct from rolling starts.
    ramp = max(.65, VMAX/(model.brake_a+float(road_load(VMAX))))
    vg = VMAX*np.clip(1.0-t/ramp,0.,1.)
    yield make('stop_wait', vg)

    # Preserve speed first and create most delay later.  The switch time is
    # chosen relative to the free-coasting arrival scale (~25 s in baseline).
    coast_until = min(.72*hi, max(1.0, min(lo if lo>0 else .42*hi, .42*hi)))
    tail = np.clip((t-coast_until)/max(.75,ramp),0.,1.)
    vg = VMAX*(1-tail) - np.minimum(t,coast_until)*float(road_load(VMAX))
    yield make('coast_then_stop', np.maximum(vg,0.))

    # Brake early, then rebuild speed around the earliest possible green. This
    # gives the narrow-prior anticipatory-acceleration basin an independent
    # initialization rather than relying on the legacy start.
    floor = min(.75*VMAX, max(1.0, .78*D0/max(hi,1e-12)))
    ramp = max(.7, (VMAX-floor)/(model.brake_a+float(road_load(VMAX))))
    vg = floor + (VMAX-floor)*np.clip(1-t/ramp,0.,1.)
    recover_start = max(0., lo - max(.4, .06*hi))
    recover_span = max(.5, hi-recover_start)
    target = min(VMAX, floor + .35*(VMAX-floor))
    vg += (target-floor)*np.clip((t-recover_start)/recover_span,0.,1.)
    yield make('brake_recover', np.minimum(vg,VMAX))

    for j, prev in enumerate(warm_starts or []):
        tv=np.asarray(prev['t'],float)
        if tv[-1] <= 0: continue
        vg=np.interp(t/hi, tv/tv[-1], np.asarray(prev['v'],float))
        mids=(t[:-1]+.5*dt)/hi
        oldm=(tv[:-1]+.5*np.diff(tv))/tv[-1]
        ug=np.interp(mids,oldm,np.asarray(prev['u'],float))
        bg=np.interp(mids,oldm,np.asarray(prev['b'],float))
        yield make(f'continuation_{j}',vg,ug,bg)


def solve_uniform(lo, hi, N=360, model=None, no_speed_increasing_throttle=False,
                  *, no_accel_before=None, multistart=True, warm_starts=None, return_candidates=True,
                  start_labels=None, max_cpu_time=12.0):
    """Best converged local candidate for a uniform switch-time prior.

    This mirrors :func:`solve_known`: several physically distinct starts plus
    optional continuation starts are optimized independently and the cheapest
    independently feasible candidate is retained.  This reduces initialization
    dependence but does not constitute a global-optimality certificate.

    ``no_speed_increasing_throttle`` forbids speed-increasing traction for
    the *entire* conditional-red branch (legacy counterfactual).
    ``no_accel_before=lo`` forbids it only up to earliest possible green.
    A restriction cutoff must coincide with a time grid node to avoid
    inadvertently constraining controls after that cutoff. In both cases,
    traction is bounded by the equilibrium cruise throttle at cell start.
    """
    if no_speed_increasing_throttle and no_accel_before is not None:
        raise ValueError("Choose full-red OR pre-green acceleration restriction, not both")
    if no_accel_before is not None:
        cutoff = float(no_accel_before)
        if not (0. <= cutoff <= lo):
            raise ValueError("Pre-green restriction must end no later than earliest green")
        if abs(cutoff * N / hi - round(cutoff * N / hi)) > 1e-8:
            raise ValueError("Pre-green cutoff must be a transcription time node")
    else:
        cutoff = None
    if not (0 <= lo < hi):
        raise ValueError('uniform prior requires 0 <= lo < hi')
    model = model or Model(300.0)
    dt = hi/N
    op = ca.Opti()
    v = op.variable(N+1); d = op.variable(N+1)
    u = op.variable(N); b = op.variable(N)
    op.subject_to(v[0] == VMAX); op.subject_to(d[0] == D0)
    op.subject_to(op.bounded(0, v, VMAX)); op.subject_to(op.bounded(0, d, D0))
    op.subject_to(op.bounded(0, u, 1)); op.subject_to(op.bounded(0, b, 1))
    obj = 0
    for k in range(N):
        if no_speed_increasing_throttle or (cutoff is not None and k*dt < cutoff-1e-9):
            op.subject_to(u[k] <= model.cruise_throttle(v[k]))
        vn, dn, stages = _rk4_step(model, v[k], d[k], u[k], b[k], dt,
                                   return_stages=True)
        for stage_v in stages:
            op.subject_to(stage_v >= 0)
        op.subject_to(v[k+1] == vn); op.subject_to(d[k+1] == dn)
        surv_dt, mass = _uniform_cell_weights(k*dt, (k+1)*dt, lo, hi)
        vm = (v[k]+v[k+1])/2; dm=(d[k]+d[k+1])/2
        obj += surv_dt*model.running_cost(vm, u[k]) + mass*model.green_value(dm, vm)
    op.minimize(obj)
    op.solver('ipopt', {'expand':True,'print_time':False},
              {'print_level':0,'max_iter':3500,'tol':1e-9,'acceptable_tol':1e-8,
               'bound_relax_factor':0,'sb':'yes','mu_strategy':'adaptive',
               'warm_start_init_point':'no','max_cpu_time':max_cpu_time})
    t=np.linspace(0,hi,N+1)
    candidates=[]; best=None
    for label,vg,dg,ug,bg in _uniform_initializations(lo,hi,N,model,warm_starts):
        if not multistart and label != 'legacy':
            break
        if start_labels is not None and label not in start_labels:
            continue
        op.set_initial(v,vg); op.set_initial(d,dg); op.set_initial(u,ug); op.set_initial(b,bg)
        try:
            sol=op.solve(); val=lambda z:np.asarray(sol.value(z)).ravel()
            V,D,U,B=val(v),val(d),val(u),val(b)
            cost=float(sol.value(obj))
            defect=0.; min_stage_speed=float(min(V))
            for k in range(N):
                vv,dd,stage=_rk4_step(model,float(V[k]),float(D[k]),
                                       float(U[k]),float(B[k]),dt,
                                       return_stages=True)
                min_stage_speed=min(min_stage_speed,*[float(z) for z in stage])
                defect=max(defect,abs(float(vv)-V[k+1]),abs(float(dd)-D[k+1]))
            feasibility=max(defect,max(0.,-min(V)),max(0.,-min(D)),
                            max(0.,max(V)-VMAX),max(0.,max(D)-D0),
                            max(0.,-min_stage_speed))
            constrained = np.arange(N)*dt < (cutoff-1e-9) if cutoff is not None else np.full(N, bool(no_speed_increasing_throttle))
            if np.any(constrained):
                feasibility=max(feasibility,float(np.max((U-np.asarray(model.cruise_throttle(V[:-1])))[constrained])))
            if not np.isfinite(cost) or feasibility > 2e-5:
                raise ValueError(f'failed independent feasibility check: {feasibility}')
            cand={'t':t,'v':V,'d':D,'u':U,'b':B,'cost':cost,
                  'T':hi,'lo':lo,'hi':hi}
            candidates.append({'start':label,'success':True,'cost':cost,
                               'min_speed':float(min(V)),
                               'min_distance':float(min(D)),
                               'min_stage_speed':min_stage_speed,
                               'max_defect':float(defect)})
            if best is None or cost < best['cost']-1e-9:
                best=cand; best['selected_start']=label
        except (RuntimeError,ValueError) as exc:
            candidates.append({'start':label,'success':False,'reason':str(exc)[:250]})
    if best is None:
        raise RuntimeError(f'No feasible uniform-prior solution U({lo},{hi}); diagnostics: {candidates}')
    if return_candidates:
        best['candidates']=candidates
    return best


def solve_uniform_family(priors, N=360, model=None, *, both_directions=True,
                         max_cpu_time=10.0, start_labels=None):
    """Continuation helper across a family of ``(lo,hi)`` uniform priors.

    Priors are ordered by midpoint for forward/reverse passes. Solutions are
    horizon-rescaled as warm starts; each pass still includes independent
    physical multistarts.
    """
    model=model or Model(300.0)
    ps=[(float(a),float(b)) for a,b in priors]
    if len(set(ps)) != len(ps):
        raise ValueError('Priors must be unique')
    key=lambda p:((p[0]+p[1])/2,p[1]-p[0])
    solved={}
    orders=(sorted(ps,key=key),sorted(ps,key=key,reverse=True)) if both_directions else (sorted(ps,key=key),)
    for order in orders:
        prev=None
        for lo,hi in order:
            warms=[prev] if prev is not None else []
            if (lo,hi) in solved: warms.append(solved[(lo,hi)])
            labels=None
            if start_labels is not None:
                labels=list(start_labels)+[f'continuation_{j}' for j in range(len(warms))]
            r=solve_uniform(lo,hi,N=N,model=model,warm_starts=warms,
                            max_cpu_time=max_cpu_time,start_labels=labels)
            if (lo,hi) not in solved or r['cost'] < solved[(lo,hi)]['cost']:
                solved[(lo,hi)]=r
            prev=solved[(lo,hi)]
    return [solved[p] for p in ps]


def integrate_red(policy, model=None, dt=0.004, until=None):
    model = model or Model(300.0)
    T = float(policy.get('T', policy['t'][-1]))
    if until is None: until=T
    n=len(policy['u']); cell=T/n
    t=0.0; v=VMAX; d=D0; ts=[t]; vs=[v]; ds=[d]
    while t < until-1e-11:
        k=min(int((t+1e-10)/cell), n-1)
        step=min(dt, until-t, (k+1)*cell-t)
        u=float(policy['u'][k]); b=float(policy['b'][k])
        def f(vv): return model.rhs_v(vv,u,b)
        k1=f(v);k2=f(v+.5*step*k1);k3=f(v+.5*step*k2);k4=f(v+step*k3)
        vnext=v+step*(k1+2*k2+2*k3+k4)/6
        dnext=d-step*(v+2*(v+.5*step*k1)+2*(v+.5*step*k2)+(v+step*k3))/6
        # Tiny direct-transcription residuals at a stopped red light should not
        # create artificial red crossing in the plotting integrator.
        if dnext < 0 and abs(dnext) < 2e-4: dnext=0.0
        if vnext < 0 and abs(vnext) < 2e-4: vnext=0.0
        t += step; v=max(0.0,vnext); d=max(0.0,dnext)
        ts.append(t);vs.append(v);ds.append(d)
    return np.array(ts),np.array(vs),np.array(ds)


def continue_green(t0, v0, d0, model=None, dt=0.004):
    model=model or Model(300.0)
    t=t0;v=v0;d=d0;ts=[t];vs=[v];ds=[d]
    while d > -POST+1e-8:
        u=1.0 if v < VMAX-1e-8 else float(model.cruise_throttle(VMAX))
        step=dt
        def f(vv): return model.rhs_v(vv,u,0.0)
        k1=f(v);k2=f(v+.5*step*k1);k3=f(v+.5*step*k2);k4=f(v+step*k3)
        vn=min(VMAX, v+step*(k1+2*k2+2*k3+k4)/6)
        dx=step*(v+vn)/2
        if d-dx < -POST:
            ratio=(d+POST)/dx; step*=ratio; vn=min(VMAX,v+step*f(v)); dx=d+POST
        t+=step;v=vn;d-=dx;ts.append(t);vs.append(v);ds.append(d)
    return np.array(ts),np.array(vs),np.array(ds)


def branch(policy, switch_time, model=None):
    model=model or Model(300.0)
    tr,vr,dr=integrate_red(policy,model,until=switch_time)
    tg,vg,dg=continue_green(switch_time,vr[-1],dr[-1],model)
    return {'t':np.r_[tr,tg[1:]], 'v':np.r_[vr,vg[1:]], 'd':np.r_[dr,dg[1:]]}


def realized_cost(policy, switch_time, model=None):
    model=model or Model(300.0)
    # Use the transcription states for the red running integral and exact green value at switch.
    t=policy['t']; u=policy['u']; v=policy['v']; d=policy['d']
    total=0.0
    for k in range(len(u)):
        a=t[k]; z=t[k+1]
        if a >= switch_time: break
        seg=min(z,switch_time)-a
        if seg <= 0: continue
        frac=seg/(z-a)
        vm=v[k] + 0.5*frac*(v[k+1]-v[k])
        total += seg*float(model.running_cost(vm,u[k]))
        if z >= switch_time: break
    vs=float(np.interp(switch_time,t,v)); ds=float(np.interp(switch_time,t,d))
    return total + float(model.green_value(ds,vs))


def save_npz(path, result):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    # Diagnostics contain dictionaries: do not introduce pickle-only object
    # arrays into canonical policy data used by the plotting scripts.
    public_keys=('t','v','d','u','b','cost','T','lo','hi','selected_start')
    np.savez(path, **{k:result[k] for k in public_keys if k in result})
