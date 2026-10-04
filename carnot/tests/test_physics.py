import math
import unittest
import numpy as np

from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.core.ccd import body_pair_toi, disc_disc_toi
from microthermo.core.boundaries import BoundaryKind, SegmentWall, Wall, resolve_wall
from microthermo.core.geometry import (closest_point_segment, convex_separation,
                                       convex_witnesses, polygon_segment_witnesses,
                                       world_polygon)
from microthermo.core.contacts import resolve_elastic
from microthermo.core.events import Contact, TOIResult, TOIStatus
from microthermo.core.mechanisms import CarnotCam, Shaft
from microthermo.core.numeric import numba_available
from microthermo.experiments.carnot import CarnotMechanism, CamPistonWall
from microthermo.core.state import BodyArrays, BodySpec
from microthermo.runner.simulation import NumericalFailure, Portal, Simulation, World
from microthermo.validation import run_validation


class PhysicsTests(unittest.TestCase):
    def test_controlled_shaft_speed_change_is_recorded_motor_intervention(self):
        config=RunConfig(preset="carnot_discs",particles=4,seed=123,
                         shaft_speed=.15)
        direct=load_preset(config)
        sampled=load_preset(config)
        for sim in (direct,sampled):
            sim.advance_to(.137)
            residual=sim.snapshot().energy_residual
            result=sim.apply_command({"name":"set_shaft_speed","speed":.3})
            self.assertTrue(result["intervention"])
            self.assertEqual(sim.events[-1].kind,"shaft_speed_change")
            self.assertAlmostEqual(sim.snapshot().energy_residual,residual,
                                   delta=1e-12)
        sampled.advance_to(.271)
        checkpoint=direct.checkpoint()
        a=direct.advance_to(.5)
        b=sampled.advance_to(.5)
        np.testing.assert_array_equal(a.position,b.position)
        self.assertEqual(direct.events,sampled.events)
        self.assertAlmostEqual(a.shaft_phase,.137*.15+(.5-.137)*.3,
                               delta=1e-12)
        direct.restore(checkpoint)
        replayed=direct.advance_to(.5)
        np.testing.assert_array_equal(replayed.position,a.position)
        self.assertEqual(replayed.shaft_phase,a.shaft_phase)
        with self.assertRaises(ValueError):
            direct.apply_command({"name":"set_shaft_speed","speed":0})

    def test_moving_thermal_wall_mode_energy_matches_heat_after_work(self):
        state=BodyArrays.from_specs([BodySpec((.9,.5),(1.,0.),radius=.1)])
        wall=Wall(np.array([1.,.5]),np.array([-1.,0.]),1.,
                  kind=BoundaryKind.HOT,temperature=1.5,
                  velocity=np.array([.1,0.]))
        diagnostics={}
        _,heat,work=resolve_wall(state,0,wall,np.array([1.,.5]),
                                  np.random.default_rng(123),
                                  contact_normal=np.array([1.,0.]),
                                  diagnostics=diagnostics)
        self.assertNotEqual(work,0.)
        self.assertAlmostEqual(
            diagnostics["outgoing_normal_mode_energy"]-
            diagnostics["incoming_normal_mode_energy"],heat,places=12)
        self.assertEqual(diagnostics["reservoir_temperature"],1.5)

    def test_scientific_validation(self):
        failed = [r for r in run_validation() if not r.passed]
        self.assertEqual([], failed)

    def test_head_on_toi(self):
        q = disc_disc_toi(np.array([0., 0.]), np.array([1., 0.]), .5,
                          np.array([3., 0.]), np.array([-1., 0.]), .5, 2.)
        self.assertEqual(q.status, TOIStatus.COLLISION)
        self.assertAlmostEqual(q.time, 1.)

    def test_unequal_mass_collision_conserves(self):
        state = BodyArrays.from_specs([BodySpec((0.,0.), (2.,.3), 2., .1),
                                       BodySpec((1.,0.), (-1.,-.4), 3., .1)])
        c = Contact(0, 1, (.5,.1), (1.,0.))
        e, p, l = state.kinetic_energy(), state.momentum(), state.angular_momentum()
        resolve_elastic(state, c)
        self.assertAlmostEqual(state.kinetic_energy(), e, places=13)
        np.testing.assert_allclose(state.momentum(), p, atol=1e-13)
        self.assertAlmostEqual(state.angular_momentum(), l, places=13)

    def test_carnot_waypoints(self):
        cam = CarnotCam.design(1., 1.3, 2., 1., 3)
        for i, area in enumerate(cam.areas):
            self.assertAlmostEqual(cam.area(i*math.pi/2), area)

    def test_free_spring_trajectory_matches_commit_and_preserves_energy(self):
        shaft=Shaft(.2,.3,2.,spring_k=18.,spring_rest=-.1)
        energy=shaft.kinetic_energy()+shaft.spring_energy()
        phi,momentum,speed=shaft.trajectory(.27)
        self.assertAlmostEqual(momentum/shaft.inertia,speed,places=14)
        shaft.advance(.27)
        self.assertAlmostEqual(shaft.phi,phi,places=14)
        self.assertAlmostEqual(shaft.momentum,momentum,places=14)
        self.assertAlmostEqual(shaft.kinetic_energy()+shaft.spring_energy(),energy,places=14)

    def test_accelerating_free_piston_collision_uses_committed_trajectory(self):
        cam=CarnotCam((1.,1.5,1.7,1.2))
        shaft=Shaft(.2,0.,1.,spring_k=100.)
        mechanism=CarnotMechanism(cam,shaft,2.,1.)
        clock=[0.]
        wall=CamPistonWall(mechanism,lambda:clock[0],1.)
        state=BodyArrays.from_specs([BodySpec((cam.piston_x(.2)-.023,.5),
                                               (0.,0.),radius=.02)])
        sim=Simulation(World(state,[wall],mechanism=mechanism,
                             metadata={"clock":clock}),max_horizon=.2)
        hit=sim._wall_toi(0,0,.2)
        self.assertEqual(hit.status,TOIStatus.COLLISION)
        self.assertGreater(hit.time,0.)
        predicted_phase=mechanism.phi_at(hit.time,0.)
        predicted_wall=wall.point_at(hit.time).copy()
        self.assertAlmostEqual(predicted_wall[0]-state.pos[0,0]-state.radius[0],
                               0.,delta=2e-10)
        snap=sim.advance_to(hit.time)
        self.assertEqual(sim.event_count,1)
        self.assertAlmostEqual(sim.world.mechanism.shaft.phi,predicted_phase,places=14)
        self.assertAlmostEqual(clock[0],hit.time,places=14)
        self.assertLess(abs(snap.energy_residual),1e-12)
        self.assertLessEqual(sim.max_penetration,sim.tol.geometry)

    def test_cam_branch_events_are_exact_and_sampling_independent(self):
        cam=CarnotCam((1.,1.5,1.7,1.2))
        def make_sim(reversed_cycle):
            direction=-1 if reversed_cycle else 1
            shaft=Shaft(direction*.2,direction,1.,prescribed_omega=direction)
            mechanism=CarnotMechanism(cam,shaft,2.,1.,reversed_cycle)
            state=BodyArrays.from_specs([BodySpec((0.,0.),(0.,0.))])
            return Simulation(World(state,[],mechanism=mechanism),max_horizon=2.)
        expected=math.pi/2-.2
        for reversed_cycle in (False,True):
            with self.subTest(reversed_cycle=reversed_cycle):
                direct=make_sim(reversed_cycle)
                direct.advance_to(1.5)
                sampled=make_sim(reversed_cycle)
                for t in (.13,.49,1.1,expected-.001,expected,1.5):
                    sampled.advance_to(t)
                self.assertEqual(len(direct.events),1)
                self.assertEqual([e.as_dict() for e in direct.events],
                                 [e.as_dict() for e in sampled.events])
                self.assertAlmostEqual(direct.events[0].time,expected,delta=1e-12)
                self.assertEqual(direct.events[0].kind,"branch_transition")
                self.assertEqual(direct.events[0].metadata["from_branch"],"hot")
                self.assertEqual(direct.events[0].metadata["to_branch"],
                                 "adiabatic_expansion")
                self.assertEqual(direct.snapshot().branch,"adiabatic_expansion")
                cycle=make_sim(reversed_cycle)
                cycle.advance_to(2*math.pi+.1)
                self.assertEqual([e.metadata["to_branch"] for e in cycle.events],
                                 ["adiabatic_expansion","cold",
                                  "adiabatic_compression","hot"])
                for index,event in enumerate(cycle.events,1):
                    self.assertAlmostEqual(event.time,index*math.pi/2-.2,
                                           delta=1e-12)

    def test_free_spring_branch_events_include_return_crossing(self):
        cam=CarnotCam((1.,1.5,1.7,1.2))
        mechanism=CarnotMechanism(cam,Shaft(.2,3.,1.,spring_k=1.),2.,1.)
        state=BodyArrays.from_specs([BodySpec((0.,0.),(0.,0.))])
        sim=Simulation(World(state,[],mechanism=mechanism),max_horizon=.2)
        sim.advance_to(3.)
        self.assertEqual([e.kind for e in sim.events],
                         ["branch_transition","branch_transition"])
        self.assertEqual([e.metadata["to_branch"] for e in sim.events],
                         ["adiabatic_expansion","hot"])
        for event in sim.events:
            self.assertAlmostEqual(event.metadata["phase"],math.pi/2,places=12)
        self.assertEqual(sim.snapshot().branch,"hot")

    def test_loaded_spring_reversals_charge_actual_angular_travel(self):
        shaft=Shaft(.4,1.,2.,spring_k=20.,load_torque=1.)
        initial=shaft.kinetic_energy()+shaft.spring_energy()
        predicted_phi,predicted_momentum=shaft.trajectory(1.5)[:2]
        turns=shaft.turning_times(1.5)
        self.assertEqual(len(turns),2)
        _,_,travel=shaft.advance_with_travel(1.5)
        self.assertAlmostEqual(shaft.phi,predicted_phi,places=14)
        self.assertAlmostEqual(shaft.momentum,predicted_momentum,places=14)
        self.assertGreater(travel,abs(shaft.phi-.4))
        self.assertAlmostEqual(initial-shaft.kinetic_energy()-shaft.spring_energy(),
                               shaft.load_torque*travel,places=13)
        split=Shaft(.4,1.,2.,spring_k=20.,load_torque=1.)
        _,_,first=split.advance_with_travel(.2)
        _,_,second=split.advance_with_travel(1.3)
        self.assertAlmostEqual(split.phi,shaft.phi,places=13)
        self.assertAlmostEqual(split.momentum,shaft.momentum,places=13)
        self.assertAlmostEqual(first+second,travel,places=13)

    def test_loaded_free_mechanism_ledger_closes_through_turns(self):
        cam=CarnotCam((1.,1.5,1.7,1.2))
        mechanism=CarnotMechanism(cam,Shaft(.4,1.,2.,spring_k=20.,
                                            load_torque=1.),2.,1.)
        state=BodyArrays.from_specs([BodySpec((0.,0.),(0.,0.))])
        sim=Simulation(World(state,[],mechanism=mechanism),max_horizon=2.)
        snap=sim.advance_to(1.5)
        self.assertLess(abs(snap.energy_residual),1e-12)
        self.assertGreater(sim.ledger.load_output.value,0.)
        self.assertAlmostEqual(sim.ledger.work_on.value,
                               -sim.ledger.load_output.value,places=13)
        self.assertEqual([e.kind for e in sim.events],
                         ["branch_transition"])

    def test_reflected_piston_inertia_force_free_flow(self):
        cam=CarnotCam((1.,1.5,1.7,1.2))
        shaft=Shaft(.4,1.,2.,piston_mass=3.,cam=cam)
        initial=shaft.kinetic_energy()
        self.assertGreater(shaft.effective_inertia(),shaft.inertia)
        predicted=shaft.trajectory(.5)
        shaft.advance(.5)
        self.assertAlmostEqual(shaft.phi,predicted[0],places=13)
        self.assertAlmostEqual(shaft.momentum,predicted[1],places=13)
        self.assertAlmostEqual(shaft.omega,predicted[2],places=13)
        self.assertAlmostEqual(shaft.kinetic_energy(),initial,places=13)
        split=Shaft(.4,1.,2.,piston_mass=3.,cam=cam)
        split.advance(.2); split.advance(.3)
        self.assertAlmostEqual(split.phi,shaft.phi,places=12)
        self.assertAlmostEqual(split.momentum,shaft.momentum,places=12)

    def test_reflected_piston_mass_closes_collision_energy(self):
        cam=CarnotCam((1.,1.5,1.7,1.2))
        shaft=Shaft(.4,0.,1.,piston_mass=2.,cam=cam)
        mechanism=CarnotMechanism(cam,shaft,2.,1.)
        clock=[0.]
        wall=CamPistonWall(mechanism,lambda:clock[0],1.)
        state=BodyArrays.from_specs([BodySpec((cam.piston_x(.4)-.12,.5),
                                               (1.,0.),radius=.02)])
        sim=Simulation(World(state,[wall],mechanism=mechanism,
                             metadata={"clock":clock}),max_horizon=.2)
        initial=sim.energy()
        hit=sim._wall_toi(0,0,.2)
        sim.advance_to(hit.time)
        self.assertEqual([e.kind for e in sim.events],["cam_recoil"])
        impulse=sim.events[0].impulse[0]
        self.assertAlmostEqual(sim.world.mechanism.shaft.momentum,
                               impulse*cam.piston_dx_dphi(.4),places=12)
        snap=sim.advance_to(.16)
        self.assertAlmostEqual(snap.energy,initial,places=12)
        self.assertLess(abs(snap.energy_residual),1e-12)
        self.assertLessEqual(sim.max_penetration,sim.tol.geometry)

    def test_controlled_reflected_mass_motor_and_load_ledger(self):
        cam=CarnotCam((1.,1.5,1.7,1.2))
        shaft=Shaft(.3,0.,2.,prescribed_omega=1.,load_torque=.2,
                    spring_k=1.,piston_mass=2.,cam=cam)
        mechanism=CarnotMechanism(cam,shaft,2.,1.)
        state=BodyArrays.from_specs([BodySpec((0.,0.),(0.,0.))])
        sim=Simulation(World(state,[],mechanism=mechanism),max_horizon=.2)
        before=sim.energy()
        snap=sim.advance_to(.5)
        self.assertAlmostEqual(sim.world.mechanism.motor_work-
                               sim.world.mechanism.load_output,
                               snap.energy-before,places=12)
        self.assertAlmostEqual(sim.ledger.work_on.value,
                               snap.energy-before,places=12)
        self.assertLess(abs(snap.energy_residual),1e-12)

    def test_controlled_piston_impact_records_motor_work(self):
        cam=CarnotCam((1.,1.5,1.7,1.2))
        mechanism=CarnotMechanism(cam,Shaft(.4,0.,1.,prescribed_omega=.2,
                                           piston_mass=2.,cam=cam),2.,1.)
        clock=[0.]
        wall=CamPistonWall(mechanism,lambda:clock[0],1.)
        state=BodyArrays.from_specs([BodySpec((cam.piston_x(.4)-.12,.5),
                                               (1.,0.),radius=.02)])
        sim=Simulation(World(state,[wall],mechanism=mechanism,
                             metadata={"clock":clock}),max_horizon=.2)
        snap=sim.advance_to(.2)
        self.assertEqual(sim.event_count,1)
        self.assertAlmostEqual(sim.world.mechanism.motor_work,
                               sim.ledger.work_on.value,places=12)
        self.assertLess(abs(snap.energy_residual),1e-12)

    def test_coupled_reflected_mass_spring_and_load_reversals(self):
        cam=CarnotCam((1.,1.5,1.7,1.2))
        shaft=Shaft(.4,1.,2.,spring_k=20.,load_torque=1.,
                    piston_mass=2.,cam=cam)
        initial=shaft.kinetic_energy()+shaft.spring_energy()
        predicted=shaft.trajectory(1.5)
        self.assertEqual(len(shaft.turning_times(1.5)),2)
        _,_,travel=shaft.advance_with_travel(1.5)
        self.assertAlmostEqual(shaft.phi,predicted[0],places=12)
        self.assertAlmostEqual(shaft.momentum,predicted[1],places=12)
        self.assertGreater(travel,abs(shaft.phi-.4))
        self.assertAlmostEqual(initial-shaft.kinetic_energy()-shaft.spring_energy(),
                               shaft.load_torque*travel,places=12)
        split=Shaft(.4,1.,2.,spring_k=20.,load_torque=1.,
                    piston_mass=2.,cam=cam)
        _,_,first=split.advance_with_travel(.2)
        _,_,second=split.advance_with_travel(1.3)
        self.assertAlmostEqual(split.kinetic_energy()+split.spring_energy()+
                               split.load_torque*(first+second),initial,places=12)
        self.assertAlmostEqual(split.phi,shaft.phi,delta=2e-3)

    def test_coupled_accelerating_piston_impact_matches_prediction(self):
        cam=CarnotCam((1.,1.5,1.7,1.2))
        mechanism=CarnotMechanism(cam,Shaft(.2,0.,1.,spring_k=100.,
                                           piston_mass=2.,cam=cam),2.,1.)
        clock=[0.]
        wall=CamPistonWall(mechanism,lambda:clock[0],1.)
        state=BodyArrays.from_specs([BodySpec((cam.piston_x(.2)-.023,.5),
                                               (0.,0.),radius=.02)])
        sim=Simulation(World(state,[wall],mechanism=mechanism,
                             metadata={"clock":clock}),max_horizon=.2)
        hit=sim._wall_toi(0,0,.2)
        self.assertEqual(hit.status,TOIStatus.COLLISION)
        predicted_phase=mechanism.phi_at(hit.time,0.)
        snap=sim.advance_to(hit.time)
        self.assertEqual(sim.event_count,1)
        self.assertAlmostEqual(sim.world.mechanism.shaft.phi,predicted_phase,places=12)
        self.assertLess(abs(snap.energy_residual),1e-12)

    def test_carnot_presets_include_reflected_piston_mass(self):
        for mode in ("controlled","free"):
            with self.subTest(mode=mode):
                sim=load_preset(RunConfig(preset="carnot_triangles",particles=16,
                                          seed=123,max_horizon=.02,shaft_mode=mode))
                self.assertGreater(sim.world.mechanism.shaft.piston_mass,0.)
                snap=sim.advance_to(.5)
                self.assertLessEqual(sim.max_penetration,sim.tol.geometry)
                self.assertLess(abs(snap.energy_residual),1e-10)

    def test_checkpoint_resume(self):
        cfg = RunConfig(preset="gas_box", particles=8, duration=.2, max_horizon=.02)
        sim = load_preset(cfg); sim.advance_to(.1); cp = sim.checkpoint()
        a = sim.advance_to(.2)
        sim.restore(cp); b = sim.advance_to(.2)
        np.testing.assert_array_equal(a.position, b.position)
        np.testing.assert_array_equal(a.velocity, b.velocity)

    def test_triangle_run_has_no_zero_time_loop(self):
        sim = load_preset(RunConfig(preset="triangle_equipartition", particles=8,
                                    duration=.3, max_horizon=.01))
        result = sim.advance_to(.3)
        self.assertEqual(result.event_count, sim.event_count)

    def test_carnot_work_ledger_closes(self):
        sim = load_preset(RunConfig(preset="carnot_discs", particles=10,
                                    duration=.5, max_horizon=.01))
        result = sim.advance_to(.5)
        self.assertLess(abs(result.energy_residual), 1e-10)

    def test_free_carnot_cam_recoil_closes(self):
        sim=load_preset(RunConfig(preset="carnot_discs",particles=10,duration=.5,
                                  max_horizon=.005,shaft_mode="free"))
        result=sim.advance_to(.5)
        self.assertLess(abs(result.energy_residual),2e-8)

    def test_scheduled_portal_energy_step(self):
        state=BodyArrays.from_specs([BodySpec((-.5,0.),(1.,0.),mass=1.,radius=.05)])
        portal=Portal(np.array([0.,0.]),np.array([-1.,0.]),1.,inside_label=1,
                      outside_label=0,delta_u=-.25,inside_energy=-.25)
        sim=Simulation(World(state,[],[portal]),seed=1,max_horizon=1.)
        before=sim.energy(); result=sim.advance_to(.75)
        self.assertEqual(result.occupancy[1],1)
        self.assertAlmostEqual(sim.energy(),before,places=13)

    def test_selective_membrane_runs(self):
        sim=load_preset(RunConfig(preset="selective_membrane",particles=10,
                                  duration=.2,max_horizon=.01))
        self.assertAlmostEqual(sim.advance_to(.2).energy_residual,0.,places=12)

    def test_triangle_wall_uses_vertex_not_circumcircle(self):
        state=BodyArrays.from_specs([BodySpec((0.,.5),(1.,0.),radius=.2,
                                               shape=1,angle=0.)])
        from microthermo.core.boundaries import Wall
        from microthermo.runner.simulation import World
        wall=Wall(np.array([1.,.5]),np.array([-1.,0.]),1.)
        sim=Simulation(World(state,[wall]),max_horizon=1.)
        q=sim._wall_toi(0,0,1.)
        self.assertEqual(q.status,TOIStatus.COLLISION)
        # A vertex points upward; the rightmost support is sqrt(3)*R/2.
        self.assertAlmostEqual(q.time,1.-math.sqrt(3.)*.1,places=8)
        self.assertAlmostEqual(q.contact.point[0],1.,places=8)

    def test_rotation_only_collision_with_disjoint_endpoints(self):
        state=BodyArrays.from_specs([
            BodySpec((0.,0.),(0.,0.),radius=.25,shape=1,omega=-5.862049951368155),
            BodySpec((.13781993863717446,-.4079521896293826),(0.,0.),
                     radius=.25,shape=1,angle=-1.9790765646962816)])
        def gap(t):
            a=world_polygon(state.polygons[0],state.pos[0],state.angle[0]+state.omega[0]*t)
            b=world_polygon(state.polygons[1],state.pos[1],state.angle[1])
            return convex_separation(a,b)[0]
        self.assertGreater(gap(0.),.01)
        self.assertGreater(gap(.2),.01)
        self.assertLess(gap(.1),0.)
        q=body_pair_toi(state,0,1,.2)
        self.assertEqual(q.status,TOIStatus.COLLISION)
        self.assertGreater(q.time,0.)
        self.assertLess(q.time,.1)
        self.assertGreaterEqual(q.contact.feature_a,0)
        self.assertGreaterEqual(q.contact.feature_b,0)

    def test_near_grazing_fixed_triangle_pair_rejects_clear_miss(self):
        # The triangles have almost parallel edges. A near miss previously
        # stalled advancement and caused thousands of horizon refinements.
        for speed in (1.,1e3):
            for offset,expected in ((.299999,TOIStatus.COLLISION),
                                    (.300001,TOIStatus.NO_COLLISION)):
                state=BodyArrays.from_specs([
                    BodySpec((-1.,0.),(speed,0.),radius=.2,shape=1),
                    BodySpec((0.,offset),(0.,0.),radius=.2,shape=1,
                             angle=1e-7)])
                for backend in (("python","numba") if numba_available()
                                else ("python",)):
                    with self.subTest(speed=speed,offset=offset,backend=backend):
                        query=body_pair_toi(state,0,1,2./speed,
                                            numeric_backend=backend)
                        self.assertEqual(query.status,expected)
                        if expected == TOIStatus.COLLISION:
                            self.assertAlmostEqual(query.time*speed,
                                                   .826794341778,places=8)
                        else:
                            simulation=Simulation(World(state.copy(),[]),
                                                  max_horizon=.1/speed,
                                                  numeric_backend=backend)
                            simulation.advance_to(2./speed)
                            self.assertEqual(simulation.event_count,0)
                            self.assertEqual(simulation.ccd_refinements,0)
                            self.assertEqual(simulation.ccd_failures,0)

    def test_polygon_witnesses_are_on_the_actual_features(self):
        a=np.array([[0.,0.],[1.,0.],[.5,1.]])
        b=np.array([[1.2,.3],[2.2,.3],[1.7,1.3]])
        lower,normal,wa,wb,fa,fb=convex_witnesses(a,b)
        self.assertGreater(lower,0.)
        self.assertAlmostEqual(np.linalg.norm(wb-wa),float((wb-wa)@normal),places=12)
        for point,polygon in ((wa,a),(wb,b)):
            distance=min(np.linalg.norm(point-closest_point_segment(point,polygon[i],
                         polygon[(i+1)%len(polygon)])[0]) for i in range(len(polygon)))
            self.assertLess(distance,1e-12)
        self.assertGreaterEqual(fa,0)
        self.assertGreaterEqual(fb,0)

    def test_high_speed_triangle_wall_contact(self):
        state=BodyArrays.from_specs([BodySpec((0.,.5),(1e5,0.),radius=.2,
                                               shape=1,angle=0.)])
        from microthermo.core.boundaries import Wall
        sim=Simulation(World(state,[Wall(np.array([1.,.5]),np.array([-1.,0.]),1.)]),
                       max_horizon=1e-4)
        q=sim._wall_toi(0,0,1e-4)
        self.assertEqual(q.status,TOIStatus.COLLISION)
        self.assertAlmostEqual(q.time,(1.-math.sqrt(3.)*.1)/1e5,places=12)

    def test_carnot_triangle_overlap_regressions(self):
        for count in (16,48):
            with self.subTest(particles=count):
                sim=load_preset(RunConfig(preset="carnot_triangles",particles=count,
                                          max_horizon=.02,seed=123))
                for k in range(20):
                    sim.advance_to((k+1)*.01)
                self.assertLessEqual(sim.max_penetration,sim.tol.geometry)
                self.assertLess(abs(sim.snapshot().energy_residual),1e-10)
                s=sim.world.bodies
                for a in range(s.n):
                    pa=world_polygon(s.polygons[a],s.pos[a],s.angle[a])
                    for b in range(a+1,s.n):
                        pb=world_polygon(s.polygons[b],s.pos[b],s.angle[b])
                        self.assertGreaterEqual(convex_separation(pa,pb)[0],-sim.tol.geometry)

    def test_carnot_triangle_sampling_does_not_change_events(self):
        config=RunConfig(preset="carnot_triangles",particles=16,seed=123,max_horizon=.02)
        direct=load_preset(config)
        direct_result=direct.advance_to(.2)
        sampled=load_preset(config)
        for time in (.007,.019,.023,.047,.061,.083,.099,.121,.147,.163,.181,.2):
            sampled_result=sampled.advance_to(time)
        self.assertEqual([(e.time,e.kind,e.participants) for e in direct.events],
                         [(e.time,e.kind,e.participants) for e in sampled.events])
        np.testing.assert_array_equal(direct_result.position,sampled_result.position)
        np.testing.assert_array_equal(direct_result.velocity,sampled_result.velocity)
        np.testing.assert_array_equal(direct_result.angle,sampled_result.angle)
        np.testing.assert_array_equal(direct_result.omega,sampled_result.omega)

    def test_checkpoint_at_intermediate_sample_keeps_pending_event(self):
        config=RunConfig(preset="carnot_triangles",particles=16,seed=123,max_horizon=.02)
        sim=load_preset(config)
        sim.advance_to(.007)
        checkpoint=sim.checkpoint()
        first=sim.advance_to(.2)
        first_events=[(e.time,e.kind,e.participants) for e in sim.events]
        sim.restore(checkpoint)
        second=sim.advance_to(.2)
        self.assertEqual(first_events,[(e.time,e.kind,e.participants) for e in sim.events])
        np.testing.assert_array_equal(first.position,second.position)
        np.testing.assert_array_equal(first.velocity,second.velocity)

    def test_indeterminate_ccd_refines_more_than_once(self):
        state=BodyArrays.from_specs([BodySpec((0.,0.),(0.,0.))])
        sim=Simulation(World(state,[]),max_horizon=.1)
        queried=[]
        def query(horizon):
            queried.append(horizon)
            return TOIResult(TOIStatus.INDETERMINATE,reason="synthetic ambiguity") if horizon>.025 else TOIResult(TOIStatus.NO_COLLISION)
        sim._earliest=query
        sim.step_collision()
        self.assertEqual(queried,[.1,.05,.025])
        self.assertEqual(sim.ccd_refinements,2)
        self.assertAlmostEqual(sim.time,.025)

    def test_simultaneous_three_disc_cluster_conserves(self):
        state=BodyArrays.from_specs([
            BodySpec((-2.,0.),(1.,0.),radius=.5),
            BodySpec((0.,0.),(0.,0.),radius=.5),
            BodySpec((2.,0.),(-1.,0.),radius=.5)])
        sim=Simulation(World(state,[]),max_horizon=2.)
        energy=sim.energy(); momentum=state.momentum().copy()
        sim.advance_to(1.1)
        self.assertEqual(sim.cluster_count,1)
        self.assertEqual(sim.event_count,2)
        np.testing.assert_allclose(sim.world.bodies.vel[:,0],[-1.,0.,1.],atol=1e-12)
        np.testing.assert_allclose(sim.world.bodies.momentum(),momentum,atol=1e-12)
        self.assertAlmostEqual(sim.energy(),energy,places=12)
        self.assertEqual(sim.events[0].time,sim.events[1].time)

    def test_near_simultaneous_cluster_is_sampling_and_checkpoint_independent(self):
        # The right-hand impact is half a time tolerance later than the left.
        # Both contacts must be resolved as one elastic cluster.
        def make_simulation():
            state=BodyArrays.from_specs([
                BodySpec((-2.,0.),(1.,0.),radius=.5),
                BodySpec((0.,0.),(0.,0.),radius=.5),
                BodySpec((2.+5e-13,0.),(-1.,0.),radius=.5)])
            return Simulation(World(state,[]),max_horizon=.2)

        direct=make_simulation()
        direct.advance_to(1.1)
        sampled=make_simulation()
        sampled.advance_to(.937)
        checkpoint=sampled.checkpoint()
        for time in (.973, .9999999999998, 1.003, 1.1):
            sampled.advance_to(time)
        first_events=[event.as_dict() for event in sampled.events]
        sampled.restore(checkpoint)
        sampled.advance_to(1.1)

        self.assertEqual(direct.cluster_count,1)
        self.assertEqual(direct.event_count,2)
        self.assertEqual(first_events,[event.as_dict() for event in sampled.events])
        self.assertEqual([event.as_dict() for event in direct.events],first_events)
        np.testing.assert_array_equal(direct.world.bodies.pos,sampled.world.bodies.pos)
        np.testing.assert_array_equal(direct.world.bodies.vel,sampled.world.bodies.vel)
        self.assertLessEqual(direct.max_penetration,direct.tol.geometry)
        self.assertAlmostEqual(direct.snapshot().energy_residual,0.,places=12)

    def test_near_simultaneous_triangle_cluster_backend_and_replay(self):
        def make_simulation(backend):
            state=BodyArrays.from_specs([
                BodySpec((-1.,0.),(1.,0.),radius=.2,shape=1,angle=.5),
                BodySpec((0.,0.),(0.,0.),radius=.2,shape=1,angle=.5),
                BodySpec((1.+5e-13,0.),(-1.,0.),radius=.2,shape=1,angle=.5)])
            return Simulation(World(state,[]),max_horizon=.1,
                              numeric_backend=backend)

        reference=make_simulation("python")
        reference.advance_to(1.)
        sampled=make_simulation("python")
        sampled.advance_to(.647)
        checkpoint=sampled.checkpoint()
        for time in (.699, .704, .831, 1.):
            sampled.advance_to(time)
        sampled_events=[event.as_dict() for event in sampled.events]
        sampled.restore(checkpoint)
        sampled.advance_to(1.)

        self.assertEqual(reference.cluster_count,1)
        self.assertEqual(reference.event_count,2)
        self.assertEqual(sampled_events,[event.as_dict() for event in sampled.events])
        self.assertEqual([event.as_dict() for event in reference.events],sampled_events)
        np.testing.assert_array_equal(reference.world.bodies.pos,sampled.world.bodies.pos)
        np.testing.assert_array_equal(reference.world.bodies.vel,sampled.world.bodies.vel)
        self.assertLessEqual(reference.max_penetration,reference.tol.geometry)
        self.assertAlmostEqual(reference.snapshot().energy_residual,0.,places=12)

        if numba_available():
            compiled=make_simulation("numba")
            compiled.advance_to(1.)
            self.assertEqual(compiled.cluster_count,reference.cluster_count)
            self.assertEqual([(e.kind,e.participants,e.time) for e in compiled.events],
                             [(e.kind,e.participants,e.time) for e in reference.events])
            np.testing.assert_allclose(compiled.world.bodies.pos,
                                       reference.world.bodies.pos,atol=1e-12)
            np.testing.assert_allclose(compiled.world.bodies.vel,
                                       reference.world.bodies.vel,atol=1e-12)
            self.assertLessEqual(compiled.max_penetration,compiled.tol.geometry)
            self.assertAlmostEqual(compiled.snapshot().energy_residual,0.,places=12)

    def test_simultaneous_corner_wall_cluster_conserves(self):
        state=BodyArrays.from_specs([BodySpec((.5,.5),(-1.,-1.),radius=.1)])
        walls=[Wall(np.array([0.,0.]),np.array([1.,0.]),1.),
               Wall(np.array([0.,0.]),np.array([0.,1.]),1.)]
        sim=Simulation(World(state,walls),max_horizon=1.)
        energy=sim.energy()
        sim.advance_to(.5)
        self.assertEqual(sim.cluster_count,1)
        np.testing.assert_allclose(sim.world.bodies.vel[0],[1.,1.],atol=1e-12)
        self.assertAlmostEqual(sim.energy(),energy,places=12)
        self.assertLessEqual(sim.max_penetration,sim.tol.geometry)

    def test_rank_deficient_duplicate_wall_cluster(self):
        state=BodyArrays.from_specs([BodySpec((.5,.5),(-1.,0.),radius=.1)])
        walls=[Wall(np.array([0.,0.]),np.array([1.,0.]),1.) for _ in range(2)]
        sim=Simulation(World(state,walls),max_horizon=1.)
        sim.advance_to(.5)
        self.assertEqual(sim.cluster_count,1)
        self.assertEqual(sim.events[0].metadata["cluster_rank"],1)
        self.assertAlmostEqual(sim.world.bodies.vel[0,0],1.,places=12)
        self.assertAlmostEqual(sim.snapshot().energy_residual,0.,places=12)

    def test_checkpoint_restores_history_after_branching(self):
        state=BodyArrays.from_specs([
            BodySpec((-2.,0.),(1.,0.),radius=.5),
            BodySpec((0.,0.),(0.,0.),radius=.5),
            BodySpec((2.,0.),(-1.,0.),radius=.5)])
        sim=Simulation(World(state,[]),max_horizon=2.)
        sim.advance_to(.5); earlier=sim.checkpoint()
        sim.advance_to(1.1); later=sim.checkpoint()
        sim.advance_to(1.5)
        original=[e.as_dict() for e in sim.events]
        sim.restore(earlier); sim.advance_to(.75)
        sim.restore(later); sim.advance_to(1.5)
        self.assertEqual(original,[e.as_dict() for e in sim.events])
        self.assertEqual(sim.cluster_count,1)

    def test_simultaneous_thermal_corner_pauses_with_checkpoint(self):
        state=BodyArrays.from_specs([BodySpec((.5,.5),(-1.,-1.),radius=.1)])
        walls=[Wall(np.array([0.,0.]),np.array([1.,0.]),1.,
                    kind=BoundaryKind.HOT,temperature=2.),
               Wall(np.array([0.,0.]),np.array([0.,1.]),1.)]
        sim=Simulation(World(state,walls),max_horizon=1.)
        with self.assertRaises(NumericalFailure):
            sim.advance_to(.5)
        self.assertEqual(sim.time,0.)
        self.assertEqual(sim.event_count,0)
        self.assertIsNotNone(sim.failure_checkpoint)
        self.assertEqual(sim.failure_diagnostic["reason"],
                         "unsupported simultaneous moving or thermal contact")

    def test_triangle_finite_segment_uses_actual_geometry(self):
        state=BodyArrays.from_specs([BodySpec((0.,0.),(1.,0.),radius=.2,shape=1)])
        wall=SegmentWall((1.,-.5),(1.,.5),thickness=.05)
        sim=Simulation(World(state,[wall]),max_horizon=1.)
        q=sim._wall_toi(0,0,1.)
        self.assertEqual(q.status,TOIStatus.COLLISION)
        expected=1.-.05-math.sqrt(3.)*.1
        self.assertAlmostEqual(q.time,expected,places=8)
        self.assertAlmostEqual(q.contact.point[0],.95,places=8)
        self.assertGreaterEqual(q.contact.feature_a,0)
        energy=sim.energy()
        sim.advance_to(1.)
        self.assertGreaterEqual(sim.event_count,1)
        self.assertAlmostEqual(sim.events[0].time,expected,places=8)
        self.assertAlmostEqual(sim.energy(),energy,places=11)
        self.assertLessEqual(sim.max_penetration,sim.tol.geometry)

    def test_triangle_misses_short_segment_outside_its_real_shape(self):
        state=BodyArrays.from_specs([BodySpec((-.5,.3),(1.,0.),radius=.2,shape=1)])
        wall=SegmentWall((0.,-.1),(0.,.1),thickness=.05)
        sim=Simulation(World(state,[wall]),max_horizon=1.)
        self.assertEqual(sim._wall_toi(0,0,1.).status,TOIStatus.NO_COLLISION)
        sim.advance_to(1.)
        self.assertEqual(sim.event_count,0)

    def test_polygon_segment_crossing_and_initial_overlap(self):
        square=np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]])
        distance,*_=polygon_segment_witnesses(square,np.array([-2.,0.]),np.array([2.,0.]))
        self.assertEqual(distance,0.)
        state=BodyArrays.from_specs([BodySpec((.95,0.),(0.,0.),radius=.2,shape=1)])
        with self.assertRaisesRegex(ValueError,"overlapping initial state"):
            Simulation(World(state,[SegmentWall((1.,-.5),(1.,.5),thickness=.05)]))

    def test_rotation_only_triangle_segment_collision_between_clear_endpoints(self):
        state=BodyArrays.from_specs([BodySpec((0.,0.),(0.,0.),radius=.25,
                                               shape=1,omega=-math.pi)])
        wall=SegmentWall((.23,-.05),(.23,.05),thickness=.005)
        for t in (0.,1.):
            polygon=world_polygon(state.polygons[0],state.pos[0],state.angle[0]+state.omega[0]*t)
            self.assertGreater(polygon_segment_witnesses(polygon,wall.start,wall.end)[0]-wall.thickness,0.)
        sim=Simulation(World(state,[wall]),max_horizon=1.)
        q=sim._wall_toi(0,0,1.)
        self.assertEqual(q.status,TOIStatus.COLLISION)
        self.assertGreater(q.time,0.)
        self.assertLess(q.time,.5)

    def test_triangle_segment_changing_closest_feature_replays(self):
        wall=SegmentWall((0.,-.2),(0.,.2),thickness=.01)
        def make_simulation():
            state=BodyArrays.from_specs([BodySpec(
                (-.559,.159),(1.128,-.255),radius=.12,shape=1,
                angle=1.725,omega=-3.035)])
            return Simulation(World(state,[wall]),max_horizon=.05)

        direct=make_simulation()
        state=direct.world.bodies
        hit=direct._wall_toi(0,0,.5)
        self.assertEqual(hit.status,TOIStatus.COLLISION)
        self.assertAlmostEqual(hit.time,.380435527208,places=9)

        def closest_feature(t):
            vertices=world_polygon(state.polygons[0],
                                   state.pos[0]+state.vel[0]*t,
                                   state.angle[0]+state.omega[0]*t)
            return polygon_segment_witnesses(vertices,wall.start,wall.end)[-2:]

        self.assertEqual(closest_feature(0.),(2,1))
        self.assertEqual(closest_feature(.8*hit.time),(4,1))
        self.assertEqual((hit.contact.feature_a,hit.contact.feature_b),(4,1))
        direct.advance_to(.5)

        sampled=make_simulation()
        sampled.advance_to(.173)
        checkpoint=sampled.checkpoint()
        for time in (.297,.369,.385,.5):
            sampled.advance_to(time)
        sampled_events=[event.as_dict() for event in sampled.events]
        sampled.restore(checkpoint)
        sampled.advance_to(.5)

        self.assertEqual(direct.event_count,1)
        self.assertEqual([event.as_dict() for event in direct.events],sampled_events)
        self.assertEqual([event.as_dict() for event in sampled.events],sampled_events)
        np.testing.assert_array_equal(direct.world.bodies.pos,sampled.world.bodies.pos)
        np.testing.assert_array_equal(direct.world.bodies.vel,sampled.world.bodies.vel)
        self.assertLessEqual(direct.max_penetration,direct.tol.geometry)
        self.assertAlmostEqual(direct.snapshot().energy_residual,0.,places=12)

    def test_high_speed_triangle_segment_contact(self):
        state=BodyArrays.from_specs([BodySpec((0.,0.),(1e5,0.),radius=.2,shape=1)])
        sim=Simulation(World(state,[SegmentWall((1.,-.5),(1.,.5),thickness=.05)]),
                       max_horizon=1e-4)
        q=sim._wall_toi(0,0,1e-4)
        self.assertEqual(q.status,TOIStatus.COLLISION)
        self.assertAlmostEqual(q.time,(1.-.05-math.sqrt(3.)*.1)/1e5,places=12)

    def test_oblique_triangle_impact_reverses(self):
        state=BodyArrays.from_specs([
            BodySpec((-.5,.08),(1.,0.),radius=.2,shape=1,angle=.2),
            BodySpec((.5,-.08),(-1.,0.),radius=.2,shape=1,angle=-.3)])
        initial=state.copy()
        sim=Simulation(World(state,[]),max_horizon=.05)
        sim.advance_to(.7)
        self.assertEqual(sim.event_count,1)
        sim.apply_command({"name":"reverse_particle_velocities"})
        sim.advance_to(1.4)
        self.assertEqual(sim.event_count,2)
        np.testing.assert_allclose(sim.world.bodies.pos,initial.pos,atol=1e-9)
        np.testing.assert_allclose(sim.world.bodies.angle,initial.angle,atol=1e-9)
        np.testing.assert_allclose(sim.world.bodies.vel,-initial.vel,atol=1e-9)
        np.testing.assert_allclose(sim.world.bodies.omega,-initial.omega,atol=1e-9)

    def test_symmetric_two_feature_triangle_impact_replays_and_reverses(self):
        def make_simulation(backend):
            state=BodyArrays.from_specs([
                BodySpec((-.5,0.),(1.,0.),radius=.2,shape=1,angle=math.pi),
                BodySpec((.5,0.),(-1.,0.),radius=.2,shape=1,angle=0.)])
            return Simulation(World(state,[]),max_horizon=.05,
                              numeric_backend=backend)

        direct=make_simulation("python")
        initial=direct.world.bodies.copy()
        direct.advance_to(.8)
        self.assertEqual(direct.event_count,4)
        self.assertEqual(direct.cluster_count,2)
        self.assertEqual(len({event.time for event in direct.events}),2)
        self.assertLessEqual(direct.max_penetration,direct.tol.geometry)
        self.assertAlmostEqual(direct.energy(),initial.kinetic_energy(),places=12)

        sampled=make_simulation("python")
        sampled.advance_to(.2)
        checkpoint=sampled.checkpoint()
        for time in (.4,.6,.8):
            sampled.advance_to(time)
        sampled.restore(checkpoint)
        sampled.advance_to(.8)
        self.assertEqual([event.as_dict() for event in direct.events],
                         [event.as_dict() for event in sampled.events])
        np.testing.assert_array_equal(direct.world.bodies.pos,sampled.world.bodies.pos)
        np.testing.assert_array_equal(direct.world.bodies.vel,sampled.world.bodies.vel)

        direct.apply_command({"name":"reverse_particle_velocities"})
        direct.advance_to(1.6)
        self.assertEqual(direct.event_count,8)
        np.testing.assert_allclose(direct.world.bodies.pos,initial.pos,atol=1e-10)
        np.testing.assert_allclose(np.sin(direct.world.bodies.angle),
                                   np.sin(initial.angle),atol=1e-10)
        np.testing.assert_allclose(np.cos(direct.world.bodies.angle),
                                   np.cos(initial.angle),atol=1e-10)
        np.testing.assert_allclose(direct.world.bodies.vel,-initial.vel,atol=1e-10)
        np.testing.assert_allclose(direct.world.bodies.omega,-initial.omega,atol=1e-10)

        if numba_available():
            compiled=make_simulation("numba")
            compiled.advance_to(.8)
            self.assertEqual([(event.kind,event.participants) for event in compiled.events],
                             [(event.kind,event.participants) for event in sampled.events])
            np.testing.assert_allclose([event.time for event in compiled.events],
                                       [event.time for event in sampled.events],atol=1e-12)
            np.testing.assert_allclose(compiled.world.bodies.pos,
                                       sampled.world.bodies.pos,atol=1e-12)
            np.testing.assert_allclose(compiled.world.bodies.vel,
                                       sampled.world.bodies.vel,atol=1e-12)
            np.testing.assert_allclose(compiled.world.bodies.omega,
                                       sampled.world.bodies.omega,atol=1e-12)


if __name__ == "__main__":
    unittest.main()
