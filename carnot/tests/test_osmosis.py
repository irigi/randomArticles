"""Osmosis presets and the first membrane experiments (osmosis plan, milestone 5)."""

import math
import unittest

import numpy as np

from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.core import edmd
from microthermo.experiments import preset_names
from microthermo.experiments.osmosis import DISC_RADIUS, design_osmosis
from microthermo.measurements.osmosis import OsmosisInstruments


def run(preset="osmosis", duration=400., settle=50., **options):
    sim = load_preset(RunConfig(preset=preset, **options))
    instruments = OsmosisInstruments(average_after=settle)
    mu = []
    for k in range(1, int(duration/.5) + 1):
        frame = instruments.observe(sim, sim.advance_to(.5*k))
        if frame.current.time > settle:
            mu.append((frame.current.mu_right, frame.current.mu_free_left))
    return sim, frame, np.nanmean(np.asarray(mu), axis=0)


class OsmosisDesignTests(unittest.TestCase):
    def test_sizing_rule_and_membrane_clearances(self):
        design = design_osmosis(RunConfig(preset="osmosis"))
        ring, membrane = design.ring, design.membrane
        self.assertAlmostEqual(8*math.pi*ring.outer_radius**2, .2*2*2)
        self.assertAlmostEqual(ring.outer_radius - ring.inner_radius, .03)
        self.assertAlmostEqual(ring.well_radius, ring.inner_radius - 2*DISC_RADIUS)
        self.assertAlmostEqual(ring.mouth_width, 6*DISC_RADIUS)
        self.assertEqual(ring.well_depth, 1.5)
        self.assertTrue(membrane.passes(DISC_RADIUS))
        self.assertFalse(membrane.passes(ring.outer_radius))
        self.assertEqual(design_osmosis(RunConfig(preset="osmosis_hard")).ring.well_depth, 0.)
        # More hosts shrink them until the well is under two disc diameters.
        self.assertLess(design_osmosis(RunConfig(preset="osmosis", hosts=10)).ring.outer_radius,
                        ring.outer_radius)
        with self.assertRaisesRegex(ValueError, "fewer hosts"):
            design_osmosis(RunConfig(preset="osmosis", hosts=16))

    def test_presets_configuration_and_engine(self):
        names = preset_names()
        self.assertIn("osmosis", names)
        self.assertIn("osmosis_hard", names)
        for old in ("labyrinth", "labyrinth_energetic", "selective_membrane"):
            self.assertNotIn(old, names)
        config = RunConfig(preset="osmosis")
        self.assertEqual((config.particles, config.resolved_engine), (200, "edmd"))
        self.assertEqual(RunConfig().resolved_engine, "reference")
        with self.assertRaises(ValueError):
            RunConfig(preset="osmosis", engine="reference").validate()
        with self.assertRaises(ValueError):
            RunConfig(preset="carnot_discs", hosts=4).validate()
        with self.assertRaises(ValueError):
            RunConfig(preset="osmosis_hard", binding_energy=1.).validate()
        # Fields added for osmosis do not change older run digests.
        self.assertEqual(RunConfig().digest(), RunConfig(engine="reference").digest())
        self.assertNotEqual(RunConfig(preset="osmosis").digest(),
                            RunConfig(preset="osmosis", hosts=6).digest())

    @unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
    def test_world_starts_with_discs_right_and_hosts_left(self):
        sim = load_preset(RunConfig(preset="osmosis"))
        world = sim.world
        hosts = sorted(world.rings)
        self.assertEqual(hosts, list(range(8)))
        x = world.bodies.pos[:, 0]
        self.assertTrue(np.all(x[:8] < world.metadata["membrane_x"]))
        self.assertTrue(np.all(x[8:] > world.metadata["membrane_x"]))
        self.assertTrue(all(w.kind.value == "hot" and w.temperature == 1. for w in world.walls))
        mixed = load_preset(RunConfig(preset="osmosis", discs_start="mixed"))
        self.assertGreater(np.sum(mixed.world.bodies.pos[8:, 0] < 2.), 20)


@unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
class MembraneExperimentTests(unittest.TestCase):
    def test_e1_free_densities_equalise_across_the_membrane(self):
        hard, hard_frame, hard_mu = run("osmosis_hard")
        bound, bound_frame, bound_mu = run("osmosis")
        for mu in (hard_mu, bound_mu):
            # mu/T on the right and of the free discs on the left agree.
            self.assertAlmostEqual(mu[0], mu[1], delta=.08)
        # Hard-only hosts follow the ideal reference closely: fewer discs on
        # the host side, because a ring takes more area than its cavity returns.
        averages, reference = hard_frame.averages, hard_frame.reference
        self.assertAlmostEqual(averages["left"], reference["left"], delta=.05*reference["left"])
        self.assertLess(averages["left"], averages["right"])
        # Binding pulls discs through the membrane onto the host side.
        self.assertGreater(bound_frame.averages["left"], averages["left"] + 4)
        for sim in (hard, bound):
            snap = sim.snapshot()
            self.assertLessEqual(abs(snap.energy_residual), 1e-10*abs(snap.energy))

    def test_e2_binding_isotherm_and_its_dilute_limit(self):
        ratios, occupancy = [], []
        for discs in (100, 400):
            sim, frame, _ = run("osmosis", 500., particles=discs)
            averages = frame.averages
            meta = sim.world.metadata
            per_host = averages["bound"]/meta["hosts"]
            c_right = averages["right"]/meta["accessible_area_right"]
            ratios.append(per_host/c_right)
            occupancy.append(per_host)
        ideal = meta["well_area"]*math.exp(1.5)
        # Crowding of finite discs lowers the ratio as the wells fill ...
        self.assertLess(ratios[1], ratios[0])
        self.assertLess(ratios[0], ideal)
        # ... and its linear extrapolation to empty wells meets the ideal value.
        slope = (ratios[1] - ratios[0])/(occupancy[1] - occupancy[0])
        dilute = ratios[0] - slope*occupancy[0]
        self.assertAlmostEqual(dilute, ideal, delta=.12*ideal)


if __name__ == "__main__":
    unittest.main()
