"""Independent trajectory and complete-cycle checks for the Carnot shaft."""

import math
import unittest

from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.core.mechanisms import CarnotCam, Shaft
from microthermo.core.state import BodyArrays, BodySpec
from microthermo.experiments.carnot import CarnotMechanism
from microthermo.runner.simulation import Simulation, World
from microthermo.measurements.cycles import summarize_cycle
from microthermo.measurements.live import LiveInstruments
from microthermo.measurements.pressure_area import compare_pressure_area


def reference_trajectory(cam, phi, momentum, duration, step):
    """Conventional RK4 on the canonical Hamilton equations, without Shaft flow."""
    inertia, piston_mass, spring_k, load = 2., 2., 20., .1

    def derivative(q, p):
        sector, u, width = cam._sector(q)
        span = cam.areas[(sector + 1) % 4] - cam.areas[sector]
        slope = cam.piston_dx_dphi(q)
        curvature = span * 60 * u * (1-u) * (1-2*u) / (width**2 * cam.height)
        mass = inertia + piston_mass * slope**2
        mass_derivative = 2 * piston_mass * slope * curvature
        return (p / mass,
                p*p * mass_derivative / (2*mass*mass) - spring_k*q - load)

    count = round(duration / step)
    h = duration / count
    for _ in range(count):
        a, b = derivative(phi, momentum)
        c, d = derivative(phi+h*a/2, momentum+h*b/2)
        e, f = derivative(phi+h*c/2, momentum+h*d/2)
        g, j = derivative(phi+h*e, momentum+h*f)
        phi += h*(a+2*c+2*e+g)/6
        momentum += h*(b+2*d+2*f+j)/6
    return phi, momentum


def integrated_shaft(cam, step, duration):
    shaft = Shaft(.4, 1., 2., spring_k=20., load_torque=.1,
                  piston_mass=2., cam=cam)
    elapsed = 0.
    while elapsed < duration-1e-14:
        interval = min(step, duration-elapsed)
        phi, momentum, _, turn = shaft._coupled_step(
            shaft.phi, shaft.momentum, interval)
        shaft.phi, shaft.momentum = phi, momentum
        elapsed += interval if turn is None else turn
    return shaft.phi, shaft.momentum


class CarnotCycleTests(unittest.TestCase):
    def test_full_gas_cycles_close_forward_and_reversed_ledgers(self):
        duration=2*math.pi/.15
        for preset,count in (("carnot_discs",2),("carnot_triangles",4)):
            for reversed_cycle in (False,True):
                with self.subTest(preset=preset,reversed_cycle=reversed_cycle):
                    sim=load_preset(RunConfig(
                        preset=preset,particles=count,seed=123,max_horizon=.1,
                        reversed_cycle=reversed_cycle))
                    snapshot=sim.advance_to(duration)
                    self.assertEqual(len(sim.cycle_markers),2)
                    cycle=summarize_cycle(*sim.cycle_markers)
                    self.assertAlmostEqual(cycle.end_time,duration,delta=2e-9)
                    self.assertAlmostEqual(cycle.end_phase,
                                           (-1 if reversed_cycle else 1)*2*math.pi,
                                           delta=2e-9)
                    self.assertGreater(cycle.event_count,20)
                    self.assertTrue(any(event.kind=="hot" for event in sim.events))
                    self.assertTrue(any(event.kind=="cold" for event in sim.events))
                    for event in (event for event in sim.events
                                  if event.kind in ("hot","cold")):
                        metadata=event.metadata
                        self.assertAlmostEqual(
                            metadata["outgoing_normal_mode_energy"]-
                            metadata["incoming_normal_mode_energy"],
                            event.heat_into_system,delta=1e-11)
                        self.assertEqual(metadata["reservoir_temperature"],
                                         sim.world.metadata[
                                             "T_hot" if event.kind=="hot" else "T_cold"])
                    self.assertTrue(any((event.metadata or {}).get("boundary")==3
                                        for event in sim.events))
                    self.assertLess(abs(cycle.total_first_law_residual),1e-10)
                    self.assertLess(abs(cycle.gas_first_law_residual),1e-10)
                    self.assertLess(abs(snapshot.energy_residual),1e-10)
                    self.assertLessEqual(sim.max_penetration,sim.tol.geometry)
                    instrument=LiveInstruments().observe(sim,snapshot).current
                    self.assertEqual(instrument.completed_cycles,1)
                    self.assertEqual(instrument.latest_cycle,cycle)
                    self.assertIsNone(instrument.efficiency_report.estimate)
                    self.assertEqual(instrument.efficiency_report.status,
                                     "awaiting_matched_cycles")
                    comparison=instrument.latest_pressure_area
                    self.assertIsNotNone(comparison)
                    self.assertAlmostEqual(comparison.event_work_by_gas,
                                           -cycle.piston_work_on_gas,places=12)
                    self.assertLessEqual(abs(comparison.difference),
                                         comparison.smoothing_error_bound+1e-10)
                    self.assertEqual(comparison.window,.25)

    def test_cycle_marker_survives_checkpoint_and_sampling(self):
        config=RunConfig(preset="carnot_discs",particles=2,seed=123,
                         max_horizon=.1)
        duration=2*math.pi/.15
        direct=load_preset(config)
        direct.advance_to(duration)
        sampled=load_preset(config)
        sampled.advance_to(duration/2)
        checkpoint=sampled.checkpoint()
        for index in range(11,21):
            sampled.advance_to(duration*index/20)
        markers=tuple(sampled.cycle_markers)
        self.assertEqual(markers,tuple(direct.cycle_markers))
        sampled.restore(checkpoint)
        sampled.advance_to(duration)
        self.assertEqual(markers,tuple(sampled.cycle_markers))
        direct_comparison=compare_pressure_area(direct,*direct.cycle_markers)
        sampled_comparison=compare_pressure_area(sampled,*sampled.cycle_markers)
        self.assertEqual(direct_comparison,sampled_comparison)

    def test_pressure_area_window_converges_toward_event_work(self):
        sim=load_preset(RunConfig(preset="carnot_discs",particles=2,
                                  seed=123,max_horizon=.1))
        sim.advance_to(2*math.pi/.15)
        wide=compare_pressure_area(sim,*sim.cycle_markers,window=.25)
        narrow=compare_pressure_area(sim,*sim.cycle_markers,window=.01)
        self.assertLess(abs(narrow.difference),abs(wide.difference)/10)
        self.assertLess(narrow.smoothing_error_bound,
                        wide.smoothing_error_bound/10)

    def test_controlled_ten_cycle_request_records_final_boundary(self):
        sim=load_preset(RunConfig(preset="carnot_discs",particles=2,
                                  seed=123,max_horizon=.1))
        duration=10*2*math.pi/.15
        sim.advance_to(duration)
        self.assertEqual(len(sim.cycle_markers),11)
        self.assertAlmostEqual(sim.cycle_markers[-1].time,duration,delta=1e-12)
        self.assertAlmostEqual(sim.cycle_markers[-1].phase,20*math.pi,
                               delta=1e-12)
        self.assertLess(abs(sim.snapshot().energy_residual),1e-10)

    def test_coupled_shaft_converges_to_independent_rk4(self):
        cam = CarnotCam((1., 1.5, 1.7, 1.2))
        # The reference interval stays on one side of the first turning point,
        # so the load is a smooth constant torque in the Hamilton equations.
        duration = .08
        fine = reference_trajectory(cam, .4, 1., duration, 1e-5)
        coarse_reference = reference_trajectory(cam, .4, 1., duration, 2e-5)
        self.assertLess(max(abs(a-b) for a,b in zip(fine,coarse_reference)), 1e-12)
        errors = []
        for step in (.04, .02, .01):
            result = integrated_shaft(cam, step, duration)
            errors.append(max(abs(a-b) for a,b in zip(result,fine)))
        self.assertGreater(errors[0]/errors[1], 3.)
        self.assertGreater(errors[1]/errors[2], 3.)
        self.assertLess(errors[-1], 1e-4)

    def test_controlled_complete_cycles_close_motor_and_load_ledger(self):
        cam = CarnotCam((1., 1.5, 1.7, 1.2))
        for direction in (1., -1.):
            with self.subTest(direction=direction):
                shaft = Shaft(0., 2*direction, 2.,
                              prescribed_omega=direction, spring_k=.3,
                              load_torque=.1, piston_mass=2., cam=cam)
                mechanism = CarnotMechanism(cam, shaft, 2., 1., direction < 0)
                body = BodyArrays.from_specs([BodySpec((0.,0.), (0.,0.))])
                sim = Simulation(World(body, [], mechanism=mechanism),
                                 max_horizon=.8)
                initial = sim.energy()
                snap = sim.advance_to(2*math.pi)
                mechanism = sim.world.mechanism
                self.assertAlmostEqual(mechanism.shaft.phi,
                                       direction*2*math.pi, delta=2e-12)
                self.assertEqual(len(sim.events), 4)
                self.assertTrue(all(event.kind == "branch_transition"
                                    for event in sim.events))
                self.assertAlmostEqual(mechanism.load_output, .1*2*math.pi,
                                       places=11)
                self.assertAlmostEqual(snap.energy-initial,
                                       sim.ledger.work_on.value, places=11)
                self.assertAlmostEqual(sim.ledger.work_on.value,
                                       mechanism.motor_work-mechanism.load_output,
                                       places=11)
                self.assertLess(abs(snap.energy_residual), 1e-10)

    def test_free_complete_cycle_converts_storage_to_load_output(self):
        cam = CarnotCam((1., 1.5, 1.7, 1.2))
        shaft = Shaft(0., 10., 2., load_torque=.1, piston_mass=2., cam=cam)
        lo, hi = 0., 2.
        for _ in range(34):
            middle = (lo+hi)/2
            if shaft.trajectory(middle)[0] >= 2*math.pi:
                hi = middle
            else:
                lo = middle
        mechanism = CarnotMechanism(cam, shaft, 2., 1.)
        body = BodyArrays.from_specs([BodySpec((0.,0.), (0.,0.))])
        sim = Simulation(World(body, [], mechanism=mechanism), max_horizon=.5)
        initial = sim.energy()
        snap = sim.advance_to(hi)
        mechanism = sim.world.mechanism
        self.assertAlmostEqual(mechanism.shaft.phi, 2*math.pi, delta=2e-7)
        self.assertEqual(len(sim.events), 4)
        self.assertAlmostEqual(mechanism.motor_work, 0., places=13)
        self.assertAlmostEqual(mechanism.load_output, .1*2*math.pi, delta=2e-8)
        self.assertAlmostEqual(initial-snap.energy,
                               mechanism.load_output, delta=1e-9)
        self.assertLess(abs(snap.energy_residual), 1e-9)


if __name__ == "__main__":
    unittest.main()
