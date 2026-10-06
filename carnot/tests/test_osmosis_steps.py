"""Temperature steps and the osmotic pressure (osmosis plan, milestone 6)."""

import unittest

import numpy as np

from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.core import edmd

STEPS = ((400., 2.), (800., .75))
WINDOWS = ((100., 400., .75), (500., 800., 2.), (900., 1200., .75))


def stepped_run(preset, seed=1):
    """Per temperature window: mean left-minus-right discs and bound discs."""
    sim = load_preset(RunConfig(preset=preset, temperature=.75, temperature_schedule=STEPS,
                                seed=seed, discs_start="mixed"))
    x_membrane = sim.world.metadata["membrane_x"]
    discs = np.asarray(sim.bh) < 0
    rows = []
    for start, end, temperature in WINDOWS:
        sim.advance_to(start)
        excess, bound = [], []
        for k in range(1, int((end - start)/.5) + 1):
            x = sim.advance_to(start + .5*k).position[discs, 0]
            excess.append(np.sum(x < x_membrane) - np.sum(x >= x_membrane))
            bound.append(sim.occupancy().sum())
        rows.append({"temperature": temperature, "excess": np.mean(excess),
                     "bound": np.mean(bound)})
    return sim, rows


def osmotic_pressure(preset, discs=100, hosts=8, temperature=1., duration=2100., settle=100.):
    """Total net x-force on the membrane per height, and its contact split."""
    sim = load_preset(RunConfig(preset=preset, particles=discs, hosts=hosts,
                                temperature=temperature, seed=1, discs_start="mixed"))
    height = sim.world.walls[0].length
    sim.advance_to(settle)
    host0, disc0 = (sim.post_impulses(k)[:, 0].sum() for k in ("hosts", "discs"))
    sim.advance_to(duration)
    host1, disc1 = (sim.post_impulses(k)[:, 0].sum() for k in ("hosts", "discs"))
    span = (duration - settle)*height
    ideal = hosts*temperature/sim.world.metadata["host_accessible_area"]
    return {"total": (host1 - host0 + disc1 - disc0)/span, "hosts": (host1 - host0)/span,
            "discs": (disc1 - disc0)/span, "ideal": ideal, "sim": sim}


@unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
class TemperatureStepTests(unittest.TestCase):
    def test_wall_temperature_command_and_schedule(self):
        sim = load_preset(RunConfig(preset="osmosis_hard", temperature_schedule=((5., 2.),)))
        sim.advance_to(4.)
        self.assertTrue(np.all(sim.wtemp == 1.))
        cp = sim.checkpoint()
        sim.advance_to(6.)
        self.assertTrue(np.all(sim.wtemp == 2.))
        self.assertEqual(sim.world.metadata["temperature"], 2.)
        changes = [e for e in sim.events if e.kind == "temperature_change"]
        self.assertEqual([(e.time, e.metadata["new_temperature"]) for e in changes], [(5., 2.)])
        after = sim.snapshot().position.copy()
        sim.restore(cp)
        self.assertTrue(np.all(sim.wtemp == 1.))
        np.testing.assert_array_equal(sim.advance_to(6.).position, after)
        result = sim.apply_command({"name": "set_wall_temperature", "temperature": .5})
        self.assertTrue(result["intervention"])
        self.assertTrue(np.all(sim.wtemp == .5))
        self.assertLessEqual(abs(sim.snapshot().energy_residual), 1e-11*sim.energy())
        with self.assertRaises(ValueError):
            RunConfig(preset="osmosis", temperature_schedule=((5., 1.), (2., 1.))).validate()
        with self.assertRaises(ValueError):
            RunConfig(preset="carnot_discs", temperature_schedule=((5., 1.),)).validate()

    def test_e3_heating_reverses_the_binding_imbalance(self):
        bound_sim, bound = stepped_run("osmosis")
        hard_sim, hard = stepped_run("osmosis_hard")
        # Binding pulls discs onto the host side when cold and releases them
        # past parity when hot, reversibly; hard-only hosts do not care.
        self.assertGreater(bound[0]["excess"], 8)
        self.assertLess(bound[1]["excess"], -1)
        self.assertGreater(bound[2]["excess"], 8)
        self.assertGreater(bound[0]["bound"], 1.8*bound[1]["bound"])
        self.assertAlmostEqual(bound[2]["bound"], bound[0]["bound"], delta=.15*bound[0]["bound"])
        for row in hard:
            self.assertLess(row["excess"], -5)
        self.assertAlmostEqual(hard[0]["excess"], hard[1]["excess"], delta=6)
        for sim in (bound_sim, hard_sim):
            snap = sim.snapshot()
            self.assertLessEqual(abs(snap.energy_residual), 1e-10*abs(snap.energy))

    def test_e4_osmotic_pressure_counts_the_hosts(self):
        hard = osmotic_pressure("osmosis_hard")
        bound = osmotic_pressure("osmosis")
        # Fewer discs and a doubled temperature: pressure / T stays.
        hot_dilute = osmotic_pressure("osmosis_hard", discs=50, temperature=2.)
        few = osmotic_pressure("osmosis_hard", hosts=2)
        base = hard["total"]
        # Host-host exclusion lifts it above the ideal-gas estimate.
        self.assertGreater(base, hard["ideal"])
        self.assertLess(base, 1.5*hard["ideal"])
        # Binding a quarter of the discs does not change it ...
        self.assertAlmostEqual(bound["total"], base, delta=.15*base)
        # ... nor does the number of discs; it scales with T.
        self.assertAlmostEqual(hot_dilute["total"]/2., base, delta=.15*base)
        # Four times the hosts: more than twice the pressure (ideal ratio 2.4).
        self.assertGreater(base, 2*few["total"])
        self.assertAlmostEqual(few["total"], few["ideal"], delta=.25*few["ideal"])
        # The split by contacting body is not the pressure: discs push hosts
        # onto the membrane, so the host part exceeds the total.
        self.assertGreater(hard["hosts"], base)
        self.assertLess(hard["discs"], 0)


if __name__ == "__main__":
    unittest.main()
