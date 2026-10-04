"""Controlled shaft speed and reproducible speed-study accounting."""

import json
import math
from pathlib import Path
import tempfile
import unittest
import numpy as np

from microthermo.api import load_preset
from microthermo.cli import main
from microthermo.config import RunConfig
from microthermo.core.boundaries import BoundaryKind
from microthermo.core.numeric import controlled_cam_disc_toi
from microthermo.measurements.speed_study import run_speed_study
from microthermo.measurements.stationarity import cycle_drift, gas_energy_drift


class SpeedStudyTests(unittest.TestCase):
    def test_compression_timing_changes_boundaries_and_preserves_ledger(self):
        fractions=(.25,.15,.25,.35)
        config=RunConfig(preset="carnot_discs",particles=16,seed=123,
                         shaft_speed=.15,cold_jacket=True,
                         cam_fractions=fractions)
        sim=load_preset(config)
        boundaries=sim.world.mechanism.cam.boundaries
        for observed,expected in zip(boundaries,
                                     (0.,.5*math.pi,.8*math.pi,
                                      1.3*math.pi,2*math.pi)):
            self.assertAlmostEqual(observed,expected)
        duration=2*math.pi/.15
        direct=load_preset(config)
        direct.advance_to(duration)
        for boundary in boundaries[1:-1]:
            self.assertTrue(any(event.kind=="branch_transition" and
                                abs(event.time-boundary/.15)<1e-10
                                for event in direct.events))
        sampled=load_preset(config)
        for boundary in boundaries[1:]:
            sampled.advance_to(boundary/.15)
        self.assertEqual(sampled.events,direct.events)
        self.assertEqual(sampled.ledger,direct.ledger)
        self.assertLess(abs(direct.snapshot().energy_residual),1e-10)
        self.assertLessEqual(direct.max_penetration,direct.tol.geometry)
        for bad in ((.25,.25,.25,.3),(.25,-.1,.25,.6),
                    (.25,.25,.5),(.25,.25,.25,math.nan)):
            with self.subTest(bad=bad),self.assertRaises(ValueError):
                RunConfig(preset="carnot_discs",cam_fractions=bad).validate()

    def test_compression_timing_compiled_cam_matches_python_reference(self):
        config=RunConfig(preset="carnot_discs",particles=8,seed=124,
                         shaft_speed=.15,cam_fractions=(.25,.15,.25,.35))
        sim=load_preset(config)
        wall=sim.world.walls[3]
        state=sim.world.bodies
        for phase in (.2,1.8,3.1,5.3):
            sim.world.mechanism.shaft.phi=phase
            piston_x=wall.point_at(sim.time)[0]
            for a in range(state.n):
                state.pos[a]=(piston_x-state.radius[a]-.04,.1+.08*a)
                state.vel[a]=(2.,0.)
            horizon=.05
            speed=wall.speed_bound(sim.time,sim.time+horizon)
            batch=controlled_cam_disc_toi(
                state,np.ones(state.n,dtype=np.bool_),sim.world.mechanism,
                speed,horizon,sim.tol)
            for a in range(state.n):
                scalar=sim._wall_toi(a,3,horizon,speed,wall.point_at(sim.time))
                self.assertEqual(int(batch[a,0]),{"no_collision":0,
                    "collision":1,"indeterminate":2}[scalar.status.value])
                if scalar.contact:
                    self.assertAlmostEqual(batch[a,1],scalar.time,delta=2e-11)

    def test_hot_jacket_switches_only_in_hot_sector_and_preserves_accounting(self):
        config=RunConfig(preset="carnot_discs",particles=16,seed=123,
                         shaft_speed=.15,cold_jacket=True,hot_jacket=True)
        sim=load_preset(config)
        hot_time=.25*math.pi/.15
        cold_time=1.25*math.pi/.15
        adiabatic_time=.75*math.pi/.15
        for wall in sim.world.walls[1:3]:
            self.assertEqual(wall.kind_at(hot_time),BoundaryKind.HOT)
            self.assertEqual(wall.kind_at(cold_time),BoundaryKind.COLD)
            self.assertEqual(wall.kind_at(adiabatic_time),BoundaryKind.SPECULAR)
            self.assertEqual(wall.temperature_at(hot_time),1.5)
        sim.advance_to(2*math.pi/.15)
        jacket_hot=[event for event in sim.events if event.kind=="hot" and
                    (event.metadata or {}).get("boundary") in (1,2)]
        self.assertGreater(len(jacket_hot),0)
        self.assertEqual(sim.ccd_failures,0)
        self.assertLess(abs(sim.snapshot().energy_residual),1e-10)
        self.assertLessEqual(sim.max_penetration,sim.tol.geometry)
        segmented=load_preset(config)
        for branch in range(1,5):
            segmented.advance_to(branch*math.pi/(2*.15))
        self.assertEqual(segmented.events,sim.events)
        self.assertEqual(segmented.ledger,sim.ledger)

    def test_optional_cold_jacket_switches_with_branch_and_closes_energy(self):
        config=RunConfig(preset="carnot_discs",particles=16,seed=123,
                         shaft_speed=.15,cold_jacket=True)
        sim=load_preset(config)
        hot_time=.25*math.pi/.15
        cold_time=1.25*math.pi/.15
        for wall in sim.world.walls[1:3]:
            self.assertEqual(wall.kind_at(hot_time),BoundaryKind.SPECULAR)
            self.assertEqual(wall.kind_at(cold_time),BoundaryKind.COLD)
            self.assertEqual(wall.temperature_at(cold_time),.75)
        sim.advance_to(2*math.pi/.15)
        jacket_contacts=[event for event in sim.events
                         if event.kind=="cold" and
                         (event.metadata or {}).get("boundary") in (1,2)]
        self.assertGreater(len(jacket_contacts),0)
        self.assertEqual(sim.ccd_failures,0)
        self.assertLess(abs(sim.snapshot().energy_residual),1e-10)
        self.assertLessEqual(sim.max_penetration,sim.tol.geometry)
        segmented=load_preset(config)
        for branch in range(1,5):
            segmented.advance_to(branch*math.pi/(2*.15))
        self.assertEqual(segmented.events,sim.events)
        self.assertEqual(segmented.ledger,sim.ledger)

    def test_shaft_speed_must_be_positive_and_finite(self):
        for speed in (0.,-.1,math.nan,math.inf):
            with self.subTest(speed=speed),self.assertRaises(ValueError):
                RunConfig(shaft_speed=speed).validate()

    def test_configured_speed_sets_exact_cycle_endpoints_in_both_directions(self):
        for speed,reversed_cycle in ((.075,False),(.3,True)):
            with self.subTest(speed=speed,reversed_cycle=reversed_cycle):
                config=RunConfig(preset="carnot_discs",particles=2,seed=123,
                                 shaft_speed=speed,reversed_cycle=reversed_cycle)
                sim=load_preset(config)
                self.assertAlmostEqual(sim.world.mechanism.shaft.omega,
                                       speed*(-1 if reversed_cycle else 1))
                duration=2*math.pi/speed
                sim.advance_to(duration)
                self.assertEqual(len(sim.cycle_markers),2)
                self.assertAlmostEqual(sim.cycle_markers[-1].time,duration)
                self.assertAlmostEqual(sim.cycle_markers[-1].phase,
                                       (-1 if reversed_cycle else 1)*2*math.pi)
                self.assertLess(abs(sim.snapshot().energy_residual),1e-10)

    def test_speed_study_reports_each_seed_and_withholds_unready_average(self):
        report=run_speed_study(preset="carnot_discs",speeds=(.3,),
                               seeds=(123,124),cycles=1,particles=2)
        group=report["speeds"][0]
        self.assertEqual(group["ready_seed_count"],0)
        self.assertIsNone(group["mean_measured_efficiency"])
        self.assertEqual([run["seed"] for run in group["runs"]],[123,124])
        self.assertTrue(all(run["completed_cycles"]==1 for run in group["runs"]))
        self.assertTrue(all(run["status"]=="awaiting_matched_cycles"
                            for run in group["runs"]))
        self.assertTrue(all(run["hot_heat_drift"]["status"]==
                            "insufficient_cycles" for run in group["runs"]))
        self.assertTrue(all(abs(run["first_law_residual"])<1e-10
                            for run in group["runs"]))
        for run in group["runs"]:
            branches=run["branch_diagnostics"]
            self.assertEqual([branch["branch"] for branch in branches],
                             ["hot","adiabatic_expansion","cold",
                              "adiabatic_compression"])
            self.assertAlmostEqual(sum(branch["heat_into_gas"]
                                       for branch in branches),
                                   run["cycles"][0]["heat_hot"]+
                                   run["cycles"][0]["heat_cold"]+
                                   run["cycles"][0]["heat_other"],places=12)
            self.assertAlmostEqual(branches[0]["heat_into_gas"],
                                   run["cycles"][0]["heat_hot"],places=12)
            for branch in branches:
                self.assertAlmostEqual(branch["normal_mode_heat_reconciliation"],
                                       0.,places=11)
            self.assertAlmostEqual(branches[-1]["gas_energy_end"]-
                                   branches[0]["gas_energy_start"],
                                   run["cycles"][0]["delta_gas_energy"],places=12)
        with self.assertRaises(ValueError):
            run_speed_study(preset="carnot_discs",speeds=(.3,.3),
                            seeds=(123,),cycles=1,particles=2)

    def test_cycle_boundary_energy_drift_is_a_screen_only(self):
        stable=gas_energy_drift([1.]*16,0)
        self.assertEqual(stable["status"],"shift_bounded")
        self.assertEqual(len(stable["early_batch_means"]),4)
        self.assertEqual(gas_energy_drift(list(range(16)),0)["status"],
                         "drift_signal")
        self.assertEqual(gas_energy_drift([1.]*15,0)["status"],
                         "insufficient_cycles")
        noisy=[0,0,2,2,0,0,2,2,1,1,3,3,1,1,3,3]
        self.assertEqual(cycle_drift(noisy,0)["status"],"inconclusive")

    def test_branch_sampling_preserves_event_chronology(self):
        config=RunConfig(preset="carnot_discs",particles=2,seed=123,
                         shaft_speed=.3)
        direct=load_preset(config)
        direct.advance_to(2*math.pi/.3)
        segmented=load_preset(config)
        for branch in range(1,5):
            segmented.advance_to(branch*math.pi/(2*.3))
        self.assertEqual([event.as_dict() for event in direct.events],
                         [event.as_dict() for event in segmented.events])
        self.assertEqual(direct.cycle_markers,segmented.cycle_markers)

    def test_run_cli_uses_requested_speed_for_cycle_duration(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/"run"
            self.assertEqual(main(["run","--preset","carnot_discs",
                                   "--particles","2","--shaft-speed","0.3",
                                   "--cold-jacket",
                                   "--cycles","1","--sample-interval","30",
                                   "--output",str(output)]),0)
            config=json.loads((output/"config.json").read_text())
            summary=json.loads((output/"summary.json").read_text())
            self.assertAlmostEqual(config["shaft_speed"],.3)
            self.assertTrue(config["cold_jacket"])
            self.assertAlmostEqual(config["duration"],2*math.pi/.3)
            self.assertEqual(summary["completed_cycles"],1)

    def test_speed_study_cli_writes_machine_readable_gate_state(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/"study.json"
            self.assertEqual(main(["speed-study","--speeds","0.3",
                                   "--seeds","123","--cycles","1",
                                   "--particles","2","--cold-jacket",
                                   "--output",str(output)]),0)
            study=json.loads(output.read_text())
            self.assertTrue(study["cold_jacket"])
            group=study["speeds"][0]
            self.assertEqual(group["runs"][0]["completed_cycles"],1)
            self.assertEqual(group["runs"][0]["status"],
                             "awaiting_matched_cycles")
            self.assertIsNone(group["mean_measured_efficiency"])


if __name__ == "__main__":
    unittest.main()
