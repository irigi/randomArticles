import unittest

from microthermo.api import load_preset
from microthermo.config import RunConfig
from microthermo.ui.scene import Camera, scene_bounds


class SceneCameraTests(unittest.TestCase):
    def test_complete_carnot_apparatus_fits_small_viewport(self):
        for preset in ("carnot_discs", "carnot_triangles"):
            sim = load_preset(RunConfig(preset=preset, particles=8, seed=123))
            bounds = scene_bounds(sim)
            for width,height in ((500,340),(1050,700)):
                camera = Camera(bounds,width,height)
                self.assertGreater(camera.scale,0.)
                for phase in (0., .8, 2.5, 4.7):
                    sim.world.mechanism.shaft.phi=phase
                    for part in sim.snapshot().apparatus:
                        for point in part.points:
                            x,y = camera.map(*point)
                            self.assertGreaterEqual(x,19.)
                            self.assertLessEqual(x,width-19.)
                            self.assertGreaterEqual(y,19.)
                            self.assertLessEqual(y,height-19.)

    def test_camera_keeps_equal_world_scales_and_resizes(self):
        sim = load_preset(RunConfig(preset="gas_box",particles=4,seed=123))
        bounds = scene_bounds(sim)
        for width,height in ((500,340),(900,400)):
            camera=Camera(bounds,width,height)
            origin=camera.map(0.,0.)
            horizontal=camera.map(1.,0.)
            vertical=camera.map(0.,1.)
            self.assertAlmostEqual(horizontal[0]-origin[0],camera.scale)
            self.assertAlmostEqual(origin[1]-vertical[1],camera.scale)


if __name__ == "__main__":
    unittest.main()
