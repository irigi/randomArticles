import math
import unittest

from microthermo.core.mechanisms import CarnotCam, Shaft
from microthermo.core.state import BodyArrays, BodySpec
from microthermo.experiments.carnot import CarnotMechanism, CamPistonWall
from microthermo.runner.simulation import Simulation, World


class ApparatusTests(unittest.TestCase):
    def make_sim(self, reversed_cycle=False):
        cam = CarnotCam((1., 1.5, 1.7, 1.2))
        direction = -1 if reversed_cycle else 1
        shaft = Shaft(.4*direction, .6*direction, 2.,
                      prescribed_omega=.3*direction, spring_k=.4,
                      piston_mass=2., load_torque=.1, cam=cam)
        mechanism = CarnotMechanism(cam, shaft, 2., 1., reversed_cycle)
        body = BodyArrays.from_specs([BodySpec((0.,0.), (0.,0.))])
        return Simulation(World(body, [], mechanism=mechanism), max_horizon=.2)

    def test_named_apparatus_snapshot_uses_collision_state_and_energy(self):
        sim = self.make_sim()
        first = sim.snapshot()
        for target in (.11, .31):
            snap = sim.advance_to(target)
            mechanism = sim.world.mechanism
            shaft, cam = mechanism.shaft, mechanism.cam
            parts = {part.name: part for part in snap.apparatus}
            self.assertEqual(len(parts), len(snap.apparatus))
            self.assertTrue({"cylinder", "piston", "piston_cam", "cam_follower",
                             "selector_cam", "selector_shoe", "thermal_wall",
                             "flywheel", "flywheel_spoke", "shaft_spring", "load"}
                            <= parts.keys())
            piston_x = cam.piston_x(shaft.phi)
            self.assertAlmostEqual(parts["piston"].points[0][0], piston_x)
            self.assertAlmostEqual(parts["cylinder"].points[2][0], piston_x)
            wall = CamPistonWall(mechanism, lambda: sim.time, cam.height)
            self.assertAlmostEqual(wall.point_at(sim.time)[0], piston_x)
            self.assertAlmostEqual(parts["piston"].velocity[0],
                                   cam.piston_dx_dphi(shaft.phi)*shaft.omega)
            self.assertEqual(parts["selector_shoe"].state, snap.branch)
            self.assertEqual(parts["thermal_wall"].state, snap.branch)
            self.assertAlmostEqual(parts["piston"].kinetic_energy+
                                   parts["flywheel"].kinetic_energy+
                                   parts["shaft_spring"].potential_energy,
                                   mechanism.energy(), places=12)
            self.assertAlmostEqual(parts["load"].work_output,
                                   sim.ledger.load_output.value, places=12)
            center, tip = parts["flywheel_spoke"].points
            self.assertAlmostEqual(math.dist(center,tip), .2*cam.height)
        self.assertEqual(first.apparatus[1].points[0][0],
                         sim.world.mechanism.cam.piston_x(.4))

    def test_cam_tracks_are_rigid_and_selector_follows_reversed_branch(self):
        sim = self.make_sim(reversed_cycle=True)
        first = {p.name:p for p in sim.snapshot().apparatus}
        initial_phi = sim.world.mechanism.shaft.phi
        later = sim.advance_to(5.)
        parts = {p.name:p for p in later.apparatus}
        phase_change = sim.world.mechanism.shaft.phi-initial_phi
        for name, rotation in (("piston_cam", -phase_change),
                               ("selector_cam", phase_change)):
            before = first[name].points
            after = parts[name].points
            width = max(sim.world.mechanism.cam.areas) / sim.world.mechanism.cam.height
            cx = (.53 if name == "piston_cam" else .8)*width
            cy = -.48*sim.world.mechanism.cam.height
            for index in (0, 19, 73, 128):
                x, y = before[index]
                expected = (cx+(x-cx)*math.cos(rotation)-
                            (y-cy)*math.sin(rotation),
                            cy+(x-cx)*math.sin(rotation)+
                            (y-cy)*math.cos(rotation))
                self.assertAlmostEqual(after[index][0], expected[0], places=12)
                self.assertAlmostEqual(after[index][1], expected[1], places=12)
        self.assertEqual(parts["selector_shoe"].state, later.branch)
        self.assertEqual(later.branch,
                         sim.world.mechanism.cam.branch(
                             sim.world.mechanism.shaft.phi, True))


if __name__ == "__main__":
    unittest.main()
