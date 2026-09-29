import unittest

from microthermo.core.mechanisms import CarnotCam, Shaft
from microthermo.core.state import BodyArrays, BodySpec
from microthermo.experiments.carnot import CarnotMechanism, CamPistonWall
from microthermo.measurements.live import LiveInstruments
from microthermo.runner.simulation import Simulation, World


class LiveInstrumentTests(unittest.TestCase):
    def test_piston_pressure_uses_raw_impulse_and_physical_window(self):
        cam=CarnotCam((1.,1.5,1.7,1.2),height=1.)
        mechanism=CarnotMechanism(cam,Shaft(.4,0.,1.,prescribed_omega=0.,
                                           piston_mass=2.,cam=cam),2.,1.)
        clock=[0.]
        wall=CamPistonWall(mechanism,lambda:clock[0],cam.height)
        particle=BodyArrays.from_specs([BodySpec(
            (cam.piston_x(.4)-.12,.5),(1.,0.),radius=.02)])
        sim=Simulation(World(particle,[wall],mechanism=mechanism,
                             metadata={"clock":clock,"T_hot":2.,"T_cold":1.}),
                       max_horizon=.05)
        instruments=LiveInstruments(history_limit=3,pressure_window=.25)
        initial=instruments.observe(sim).current
        self.assertIsNone(initial.pressure)
        snap=sim.advance_to(.15)
        current=instruments.observe(sim,snap).current
        self.assertEqual(len(sim.events),1)
        self.assertAlmostEqual(current.pressure,
                               abs(sim.events[0].impulse[0])/(.25*cam.height),
                               places=12)
        self.assertAlmostEqual(current.pressure_window,.25)
        self.assertAlmostEqual(current.area,cam.area(mechanism.shaft.phi))
        self.assertAlmostEqual(current.gas_energy+current.piston_energy+
                               current.flywheel_energy+current.spring_energy,
                               snap.energy,places=12)
        self.assertLess(abs(current.first_law_residual),1e-12)
        for target in (.2,.3,.4,.5):
            sim.advance_to(target)
            frame=instruments.observe(sim)
        self.assertEqual(len(frame.history),3)
        self.assertEqual(frame.history[-1].time,.5)
        self.assertAlmostEqual(frame.current.pressure_window,.25)
        self.assertAlmostEqual(frame.current.pressure,0.,places=12)


if __name__ == "__main__":
    unittest.main()
