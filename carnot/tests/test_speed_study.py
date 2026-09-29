"""Controlled shaft speed and reproducible speed-study accounting."""

import json
import math
from pathlib import Path
import tempfile
import unittest

from microthermo.api import load_preset
from microthermo.cli import main
from microthermo.config import RunConfig
from microthermo.measurements.speed_study import run_speed_study
from microthermo.measurements.stationarity import cycle_drift, gas_energy_drift


class SpeedStudyTests(unittest.TestCase):
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
                                   "--cycles","1","--sample-interval","30",
                                   "--output",str(output)]),0)
            config=json.loads((output/"config.json").read_text())
            summary=json.loads((output/"summary.json").read_text())
            self.assertAlmostEqual(config["shaft_speed"],.3)
            self.assertAlmostEqual(config["duration"],2*math.pi/.3)
            self.assertEqual(summary["completed_cycles"],1)

    def test_speed_study_cli_writes_machine_readable_gate_state(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/"study.json"
            self.assertEqual(main(["speed-study","--speeds","0.3",
                                   "--seeds","123","--cycles","1",
                                   "--particles","2","--output",str(output)]),0)
            group=json.loads(output.read_text())["speeds"][0]
            self.assertEqual(group["runs"][0]["completed_cycles"],1)
            self.assertEqual(group["runs"][0]["status"],
                             "awaiting_matched_cycles")
            self.assertIsNone(group["mean_measured_efficiency"])


if __name__ == "__main__":
    unittest.main()
