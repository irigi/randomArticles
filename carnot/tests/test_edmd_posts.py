"""Fixed circular posts and the porous membrane (osmosis plan, milestone 2)."""

import math
import unittest

import numpy as np

from microthermo.core import edmd
from microthermo.core.boundaries import BoundaryKind, CircularPost
from microthermo.core.state import BodyArrays, BodySpec
from microthermo.experiments.common import box_walls
from microthermo.experiments.membrane import build_membrane
from microthermo.runner.edmd_simulation import EdmdSimulation, edmd_unsupported_reason
from microthermo.runner.simulation import Simulation, World


def world_of(specs, posts, width=2., height=1., thermal=None):
    return World(BodyArrays.from_specs(specs), box_walls(width, height, thermal, 1.),
                 posts=list(posts))


def scatter(rng, count, radius, x0, x1, y0, y1, posts, others=(), mass=1.):
    """Random non-overlapping discs clear of the posts."""
    specs = list(others)
    placed = []
    while len(placed) < count:
        p = np.array([rng.uniform(x0 + radius, x1 - radius), rng.uniform(y0 + radius, y1 - radius)])
        if any(np.hypot(*(p - q.position)) <= radius + q.radius + 1e-3 for q in specs):
            continue
        if any(np.hypot(*(p - k.center)) <= radius + k.radius + 1e-3 for k in posts):
            continue
        spec = BodySpec(tuple(p), tuple(rng.normal(0, math.sqrt(1/mass), 2)), mass=mass,
                        radius=radius)
        specs.append(spec)
        placed.append(spec)
    return placed


def single_cell(sim):
    """Rebuild a simulation with one grid cell (every pair and post tested)."""
    x0, y0 = sim.fp[edmd.FP_LOX], sim.fp[edmd.FP_LOY]
    sim.fp[edmd.FP_CX] *= sim.ip[0]
    sim.fp[edmd.FP_CY] *= sim.ip[1]
    sim.ip[:] = 1
    sim.head = np.full(1, -1, np.int64)
    sim.cps, sim.cpi = edmd.post_cell_lists(sim.ppos, sim.prad, 1., sim.fp, sim.ip)
    edmd.build_grid(sim._bodies(), sim._grid(), sim.time)
    sim._rebuild()
    return sim


@unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
class PostAnalyticTests(unittest.TestCase):
    def test_head_on_post_time_reflection_and_impulse(self):
        post = CircularPost((1.5, .5), .1)
        sim = EdmdSimulation(world_of([BodySpec((.5, .5), (2., 0.), radius=.05)], [post]))
        sim.step_collision()
        self.assertAlmostEqual(sim.time, (1.5 - .15 - .5)/2., delta=1e-15)
        event = sim.events[-1]
        self.assertEqual((event.kind, event.metadata), ("post", {"post": 0}))
        np.testing.assert_allclose(sim.world.bodies.vel[0], [-2., 0.])
        np.testing.assert_allclose(sim.post_impulses(), [[4., 0.]])

    def test_glancing_post_reflection_follows_the_normal(self):
        post = CircularPost((1.5, .5), .1)
        b = .09            # impact parameter, less than r_disc + r_post = .15
        sim = EdmdSimulation(world_of([BodySpec((.5, .5 + b), (1., 0.), radius=.05)], [post]))
        sim.step_collision()
        p = sim.world.bodies.pos[0]
        n = (post.center - p)/np.linalg.norm(post.center - p)
        self.assertAlmostEqual(np.linalg.norm(post.center - p), .15, delta=1e-14)
        expected = np.array([1., 0.]) - 2*np.dot([1., 0.], n)*n
        np.testing.assert_allclose(sim.world.bodies.vel[0], expected, atol=1e-15)
        self.assertAlmostEqual(np.linalg.norm(sim.world.bodies.vel[0]), 1., delta=1e-15)

    def test_reference_engine_and_out_of_box_posts_are_rejected(self):
        disc = [BodySpec((.5, .5), (1., 0.), radius=.05)]
        with self.assertRaises(ValueError):
            Simulation(world_of(disc, [CircularPost((1.5, .5), .1)]))
        self.assertIn("outside", edmd_unsupported_reason(
            world_of(disc, [CircularPost((3., .5), .1)])))


@unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
class PostDynamicsTests(unittest.TestCase):
    def test_large_posts_spanning_many_cells_never_overlap(self):
        rng = np.random.default_rng(11)
        posts = [CircularPost((.6, .5), .3), CircularPost((1.45, .3), .2),
                 CircularPost((1.6, .85), .07)]
        specs = scatter(rng, 120, .015, 0., 2., 0., 1., posts)
        sim = EdmdSimulation(world_of(specs, posts))
        self.assertGreater(len(sim.cpi), 3*20)      # each post listed in many cells
        e0, p0 = sim.energy(), sim.world.bodies.momentum().copy()
        for k in range(1, 41):
            sim.advance_to(k*.25)      # every sample checks all gaps
        self.assertGreater(sum(e.kind == "post" for e in sim.events), 1000)
        self.assertLessEqual(abs(sim.energy() - e0), 1e-12*e0)
        self.assertLessEqual(sim.max_penetration, 1e-12)
        # Momentum changes only through walls and posts.
        ledger = sim.ledger
        support = np.array([ledger.support_impulse_x.value, ledger.support_impulse_y.value])
        np.testing.assert_allclose(sim.world.bodies.momentum() - p0, support, atol=1e-10)
        walls = sum(np.asarray(e.impulse) for e in sim.events
                    if e.kind == "specular")
        np.testing.assert_allclose(walls + sim.post_impulses().sum(axis=0), -support,
                                   atol=1e-10)

    def test_post_lists_agree_with_a_single_cell(self):
        rng = np.random.default_rng(5)
        posts = [CircularPost((.7, .45), .25), CircularPost((1.5, .6), .12)]
        specs = scatter(rng, 80, .02, 0., 2., 0., 1., posts)
        grid = EdmdSimulation(world_of(specs, posts))
        one = single_cell(EdmdSimulation(world_of(specs, posts)))
        self.assertGreater(grid.ip[0]*grid.ip[1], 100)
        grid.advance_to(1.)
        one.advance_to(1.)
        pairs = list(zip(grid.events, one.events))[:150]
        self.assertEqual(len(pairs), 150)
        self.assertTrue(any(a.kind == "post" for a, _ in pairs))
        for a, b in pairs:
            self.assertEqual((a.kind, a.participants, a.metadata),
                             (b.kind, b.participants, b.metadata))
            self.assertAlmostEqual(b.time, a.time, delta=1e-9)


@unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
class MembraneTests(unittest.TestCase):
    def test_builder_gives_even_gaps_and_checks_clearances(self):
        m = build_membrane(2., 0., 2., .05, .1, .02, blocked_radius=.15, mouth_width=.08)
        ys = [p.center[1] for p in m.posts]
        gaps = np.diff(ys) - 2*m.post_radius
        np.testing.assert_allclose(gaps, m.gap)
        self.assertAlmostEqual(ys[0] - m.post_radius, m.gap)
        self.assertAlmostEqual(2. - ys[-1] - m.post_radius, m.gap)
        self.assertTrue(m.passes(.02))
        self.assertFalse(m.passes(.15))
        with self.assertRaises(ValueError):          # discs cannot pass
            build_membrane(2., 0., 2., .05, .03, .02)
        with self.assertRaises(ValueError):          # hosts would pass
            build_membrane(2., 0., 2., .05, .1, .02, blocked_radius=.04)
        with self.assertRaises(ValueError):          # a post fits in a mouth
            build_membrane(2., 0., 2., .03, .1, .02, mouth_width=.08)

    def test_small_discs_cross_and_large_discs_stay(self):
        rng = np.random.default_rng(21)
        height, width = 1., 2.
        m = build_membrane(1., 0., height, .03, .08, .015, blocked_radius=.06)
        big = scatter(rng, 6, .06, 0., m.x - m.post_radius, 0., height, m.posts, mass=20.)
        small = scatter(rng, 150, .015, 0., m.x - m.post_radius, 0., height, m.posts,
                        others=big)
        sim = EdmdSimulation(world_of(big + small, m.posts, width, height,
                                      thermal=BoundaryKind.HOT))
        self.assertTrue(np.all(sim.world.bodies.pos[:, 0] < m.x))   # all start left
        right_fraction, big_hits = [], 0
        for k in range(1, 401):
            snap = sim.advance_to(k*.5)
            x = snap.position[:, 0]
            self.assertTrue(np.all(x[:6] < m.x), "a blocked disc crossed the membrane")
            right_fraction.append(np.mean(x[6:] > m.x))
            big_hits += sum(1 for e in sim.drain_events()
                            if e.kind == "post" and e.participants[0] < 6)
        self.assertGreater(big_hits, 10)           # the big discs really met the posts
        # Ends near even (the right side has a bit more free area, since the
        # big discs stay left).
        late = np.mean(right_fraction[200:])
        self.assertGreater(late, .45)
        self.assertLess(late, .6)


if __name__ == "__main__":
    unittest.main()
