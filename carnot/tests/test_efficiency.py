"""Completed-cycle efficiency gating and paired block uncertainty."""

import math
import json
import tempfile
import unittest

from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.io.exports import export_run
from microthermo.measurements.cycles import CycleSummary, assess_efficiency


def cycle(index, hot, output, direction=1, start_phase=None):
    phase=direction*index*2*math.pi if start_phase is None else start_phase
    return CycleSummary(
        start_time=float(index),end_time=float(index+1),
        start_phase=phase,end_phase=phase+direction*2*math.pi,
        heat_hot=hot,heat_cold=-hot+output,heat_other=0.,
        piston_work_on_gas=0.,motor_work=-output,load_output=0.,
        delta_gas_energy=0.,delta_apparatus_energy=0.,
        total_first_law_residual=0.,gas_first_law_residual=0.,
        event_count=10)


class EfficiencyTests(unittest.TestCase):
    def test_transients_are_excluded_and_ratio_uses_paired_cycle_sums(self):
        data=[cycle(0,100.,-50.),cycle(1,100.,-50.)]
        hot=[1.,2.,1.,2.,1.,2.,1.,2.]
        output=[.2,.7,.2,.7,.2,.7,.2,.7]
        data.extend(cycle(index+2,q,w) for index,(q,w) in
                    enumerate(zip(hot,output)))
        report=assess_efficiency(data,transient_cycles=2,min_cycles=8)
        self.assertEqual(report.status,"ready")
        estimate=report.estimate
        self.assertEqual(estimate.sample_count,8)
        self.assertEqual(estimate.block_count,4)
        self.assertEqual(estimate.direction,1)
        self.assertAlmostEqual(estimate.hot_heat,12.)
        self.assertAlmostEqual(estimate.external_output,3.6)
        self.assertAlmostEqual(estimate.value,.3)
        self.assertGreaterEqual(estimate.block_standard_error,0.)
        self.assertEqual(estimate.storage_change,0.)

    def test_waits_for_sample_count_and_positive_hot_heat_per_block(self):
        data=[cycle(i,1.,.25) for i in range(9)]
        waiting=assess_efficiency(data,transient_cycles=2,min_cycles=8)
        self.assertEqual(waiting.status,"awaiting_matched_cycles")
        self.assertIsNone(waiting.estimate)
        data.append(cycle(9,1.,.25))
        self.assertEqual(assess_efficiency(data).status,"ready")
        data[-1]=cycle(9,-10.,.25)
        rejected=assess_efficiency(data)
        self.assertEqual(rejected.status,"insufficient_hot_heat")
        self.assertIsNone(rejected.estimate)

    def test_direction_change_starts_a_new_matched_run(self):
        forward=[cycle(i,1.,.25) for i in range(10)]
        reverse=[cycle(10+i,1.,.2,direction=-1,
                       start_phase=forward[-1].end_phase-i*2*math.pi)
                 for i in range(10)]
        report=assess_efficiency(forward+reverse)
        self.assertEqual(report.status,"ready")
        self.assertEqual(report.matched_cycles,8)
        self.assertEqual(report.estimate.direction,-1)
        self.assertAlmostEqual(report.estimate.value,.2)

    def test_incomplete_run_export_withholds_efficiency(self):
        config=RunConfig(preset="carnot_discs",particles=2,seed=123,
                         max_horizon=.1)
        sim=load_preset(config)
        sim.advance_to(2*math.pi/.15)
        with tempfile.TemporaryDirectory() as directory:
            path=export_run(directory,config,sim)
            summary=json.loads((path/"summary.json").read_text())
        self.assertEqual(summary["completed_cycles"],1)
        self.assertIsNone(summary["efficiency"])
        self.assertEqual(summary["efficiency_report"]["status"],
                         "awaiting_matched_cycles")

    def test_gate_parameters_are_validated(self):
        with self.assertRaises(ValueError):
            RunConfig(transient_cycles=-1).validate()
        with self.assertRaises(ValueError):
            RunConfig(efficiency_min_cycles=3).validate()


if __name__ == "__main__":
    unittest.main()
