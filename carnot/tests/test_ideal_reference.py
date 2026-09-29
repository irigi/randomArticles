"""Analytical Carnot overlay checks independent of measured trajectories."""

import math
import unittest

from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.measurements.ideal import ideal_carnot_reference


class IdealReferenceTests(unittest.TestCase):
    def test_disc_and_triangle_branches_join_and_have_carnot_work(self):
        for preset,dof in (("carnot_discs",2),("carnot_triangles",3)):
            with self.subTest(preset=preset):
                sim=load_preset(RunConfig(preset=preset,particles=4,seed=123))
                reference=ideal_carnot_reference(sim,points_per_branch=257)
                self.assertIsNotNone(reference)
                self.assertEqual(reference.degrees_of_freedom,dof)
                self.assertEqual(reference.particles,4)
                self.assertAlmostEqual(reference.efficiency,
                    1-reference.cold_temperature/reference.hot_temperature)
                for index,branch in enumerate(reference.branches):
                    following=reference.branches[(index+1)%4]
                    self.assertAlmostEqual(branch.area[-1],following.area[0])
                    self.assertAlmostEqual(branch.pressure[-1],following.pressure[0])
                    self.assertEqual(len(branch.area),257)
                area1,area2=sim.world.mechanism.cam.areas[:2]
                self.assertAlmostEqual(reference.branches[0].pressure[0],
                                       4*reference.hot_temperature/area1)
                expected=4*(reference.hot_temperature-reference.cold_temperature)*math.log(area2/area1)
                self.assertAlmostEqual(reference.work_by_gas,expected)
                integrated=sum(sum(.5*(branch.pressure[j]+branch.pressure[j+1])*
                    (branch.area[j+1]-branch.area[j]) for j in range(256))
                    for branch in reference.branches)
                self.assertAlmostEqual(integrated,expected,delta=2e-4)

    def test_reference_does_not_depend_on_measured_state_or_cycle_direction(self):
        config=RunConfig(preset="carnot_discs",particles=4,seed=123)
        sim=load_preset(config)
        initial=ideal_carnot_reference(sim)
        sim.advance_to(.2)
        self.assertEqual(ideal_carnot_reference(sim),initial)
        reversed_sim=load_preset(RunConfig(preset="carnot_discs",particles=4,
                                            seed=123,reversed_cycle=True))
        self.assertEqual(ideal_carnot_reference(reversed_sim),initial)

    def test_non_carnot_preset_has_no_reference(self):
        sim=load_preset(RunConfig(preset="gas_box",particles=4,seed=123))
        self.assertIsNone(ideal_carnot_reference(sim))


if __name__ == "__main__":
    unittest.main()
