"""Rotating rough ring hosts in the edmd kernel (osmosis plan, milestone 3)."""

import math
import unittest

import numpy as np

from microthermo.core import edmd
from microthermo.core.boundaries import BoundaryKind, CircularPost
from microthermo.core.state import BodyArrays, BodySpec, RingGeometry, Shape
from microthermo.experiments.common import box_walls
from microthermo.runner.edmd_simulation import EdmdSimulation, edmd_unsupported_reason
from microthermo.runner.simulation import Simulation, World

RING = RingGeometry(.12, .15, .08)        # mouth: two disc diameters for r = .02


def disc(x, y, vx=0., vy=0., omega=0., r=.02, m=1.):
    return BodySpec((x, y), (vx, vy), mass=m, radius=r, omega=omega, inertia=.5*m*r*r)


def host(x, y, vx=0., vy=0., angle=0., omega=0., m=25., ring=RING):
    return BodySpec((x, y), (vx, vy), mass=m, radius=ring.outer_radius, angle=angle,
                    omega=omega, inertia=ring.uniform_inertia(m), shape=Shape.HOST)


def make(specs, rings, width=2., height=1., thermal=None, roughness=1., posts=()):
    world = World(BodyArrays.from_specs(specs), box_walls(width, height, thermal, 1.),
                  posts=list(posts), rings=rings, contact_roughness=roughness)
    return EdmdSimulation(world)


def invariants(sim):
    s = sim.world.bodies
    return sim.energy(), s.momentum().copy(), s.angular_momentum()


def gas_with_hosts(rng, hosts, discs, width, height, temperature=1., ring=RING,
                   disc_radius=.02, host_mass=25.):
    """Hosts at random angles and spins; discs outside the hosts' outer circles."""
    specs, centres = [], []
    while len(centres) < hosts:
        p = rng.uniform(ring.outer_radius + .01, [width - ring.outer_radius - .01,
                                                  height - ring.outer_radius - .01])
        if all(np.hypot(*(p - q)) > 2*ring.outer_radius + .02 for q in centres):
            centres.append(p)
            inertia = ring.uniform_inertia(host_mass)
            specs.append(host(*p, *rng.normal(0, math.sqrt(temperature/host_mass), 2),
                              angle=rng.uniform(-math.pi, math.pi),
                              omega=rng.normal(0, math.sqrt(temperature/inertia)),
                              m=host_mass, ring=ring))
    inertia = .5*disc_radius**2
    while len(specs) < hosts + discs:
        p = rng.uniform(disc_radius + .005, [width - disc_radius - .005,
                                             height - disc_radius - .005])
        if any(np.hypot(*(p - q)) < ring.outer_radius + disc_radius + .005 for q in centres):
            continue
        if any(np.hypot(*(p - s.position)) < 2*disc_radius + .005 for s in specs[hosts:]):
            continue
        specs.append(disc(*p, *rng.normal(0, math.sqrt(temperature), 2),
                          omega=rng.normal(0, math.sqrt(temperature/inertia)), r=disc_radius))
    return specs, {i: ring for i in range(hosts)}


@unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
class RingGeometryTests(unittest.TestCase):
    def test_ring_geometry_and_rejections(self):
        self.assertAlmostEqual(RING.mid, .135)
        self.assertAlmostEqual(RING.cap, .015)
        # Clear chord between the rounded ends equals the mouth width.
        beta = RING.mouth_half_angle
        self.assertAlmostEqual(2*RING.mid*math.sin(beta) - 2*RING.cap, RING.mouth_width)
        with self.assertRaises(ValueError):
            RingGeometry(.15, .12, .08)
        specs = [host(1., .5), disc(.3, .5)]
        spinless = [host(1., .5), BodySpec((.3, .5), (0., 0.), radius=.02)]
        world = World(BodyArrays.from_specs(spinless), box_walls(2., 1.), rings={0: RING})
        self.assertIn("spin", edmd_unsupported_reason(world))
        world.contact_roughness = 0.
        self.assertIsNone(edmd_unsupported_reason(world))
        with self.assertRaises(ValueError):
            Simulation(World(BodyArrays.from_specs(specs), box_walls(2., 1.), rings={0: RING}))


@unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
class RingContactTests(unittest.TestCase):
    def assert_conserved(self, before, after, tol=1e-13):
        (e0, p0, l0), (e1, p1, l1) = before, after
        self.assertAlmostEqual(e1, e0, delta=tol*max(1., e0))
        np.testing.assert_allclose(p1, p0, atol=tol)
        self.assertAlmostEqual(l1, l0, delta=tol)

    def test_radial_hit_on_outer_arc(self):
        # Mouth faces +x (angle 0); the disc comes from -x, onto the arc.
        sim = make([host(1., .5), disc(.5, .5, 1.)], {0: RING})
        before = invariants(sim)
        sim.step_collision()
        event = sim.events[-1]
        self.assertEqual((event.kind, event.participants), ("host", (1, 0)))
        self.assertEqual(event.metadata, {"feature": "outer", "rough": True})
        self.assertAlmostEqual(sim.time, (1. - .15 - .02 - .5)/1., delta=1e-15)
        # Radial: rough and smooth agree; 1-D elastic exchange with M = 25.
        v = sim.world.bodies.vel
        self.assertAlmostEqual(v[1, 0], (1 - 25)/26, delta=1e-14)
        self.assertAlmostEqual(v[0, 0], 2/26, delta=1e-14)
        self.assertEqual(sim.world.bodies.omega[0], 0.)
        self.assert_conserved(before, invariants(sim))

    def test_glancing_rough_hit_reverses_the_contact_velocity(self):
        sim = make([host(1., .5), disc(.5, .62, 1., 0., omega=3.)], {0: RING})
        before = invariants(sim)
        s = sim.world.bodies
        sim.step_collision()
        self.assertEqual(sim.events[-1].metadata["feature"], "outer")
        after = invariants(sim)
        self.assert_conserved(before, after)
        self.assertNotEqual(s.omega[0], 0.)            # the host now spins
        # Relative velocity of the contact points before and after.
        p, h = s.pos[1], s.pos[0]
        n = (h - p)/np.linalg.norm(h - p)
        c = p + .02*n

        def contact_velocity(vd, wd, vh, wh):
            rd, rh = c - p, c - h
            ud = vd + wd*np.array([-rd[1], rd[0]])
            uh = vh + wh*np.array([-rh[1], rh[0]])
            return uh - ud
        incoming = contact_velocity(np.array([1., 0.]), 3., np.zeros(2), 0.)
        outgoing = contact_velocity(s.vel[1], s.omega[1], s.vel[0], s.omega[0])
        np.testing.assert_allclose(outgoing, -incoming, atol=1e-13)

    def test_smooth_contacts_never_spin_the_host(self):
        sim = make([host(1., .5), disc(.5, .62, 1., 0., omega=3.)], {0: RING}, roughness=0.)
        sim.step_collision()
        # Arc normals pass through the centre: no torque beyond rounding.
        self.assertAlmostEqual(sim.world.bodies.omega[0], 0., delta=1e-15)
        self.assertEqual(sim.world.bodies.omega[1], 3.)
        self.assertEqual(sim.events[-1].metadata["rough"], False)

    def test_disc_enters_through_the_mouth_and_hits_the_inner_arc(self):
        # Mouth faces -x toward the disc (angle pi); aim at the centre.
        sim = make([host(1., .5, angle=math.pi), disc(.5, .5, 1.)], {0: RING})
        sim.step_collision()
        event = sim.events[-1]
        self.assertEqual(event.metadata["feature"], "inner")
        # Through the centre to the far inner surface at x = 1 + .12 - .02.
        self.assertAlmostEqual(sim.time, (1. + .12 - .02) - .5, delta=1e-13)
        x = sim.world.bodies.pos[1]
        self.assertAlmostEqual(np.linalg.norm(x - sim.world.bodies.pos[0]), .10, delta=1e-12)

    def test_rotating_cap_sweeps_into_a_disc_in_the_mouth(self):
        # Disc parked in the mouth band on the axis; the host turns, so the
        # rounded end at +beta sweeps onto it (minus-cap leads for omega < 0).
        sim = make([host(1., .5, angle=0., omega=-2.), disc(1.135, .5, 0., 0.)], {0: RING})
        before = invariants(sim)
        sim.step_collision()
        event = sim.events[-1]
        self.assertEqual(event.kind, "host")
        self.assertIn(event.metadata["feature"], ("cap_plus", "cap_minus"))
        s = sim.world.bodies
        beta = RING.mouth_half_angle
        # Contact when the cap centre is .015 + .02 from the disc centre.
        sign = 1. if event.metadata["feature"] == "cap_plus" else -1.
        cap_angle = s.angle[0] + sign*beta
        cap = s.pos[0] + RING.mid*np.array([math.cos(cap_angle), math.sin(cap_angle)])
        self.assertAlmostEqual(np.linalg.norm(s.pos[1] - cap), .035, delta=1e-9)
        self.assert_conserved(before, invariants(sim), tol=1e-12)
        self.assertGreater(np.linalg.norm(s.vel[1]), 0.)

    def test_host_host_wall_and_post_contacts(self):
        sim = make([host(.5, .5, 1., 0., omega=1.), host(1.2, .52, -1., 0.)],
                   {0: RING, 1: RING})
        before = invariants(sim)
        sim.step_collision()
        self.assertEqual(sim.events[-1].kind, "host_host")
        self.assertAlmostEqual(np.linalg.norm(sim.world.bodies.pos[1] - sim.world.bodies.pos[0]),
                               .30, delta=1e-13)
        self.assert_conserved(before, invariants(sim))
        # Wall: normal impulse through the centre leaves the spin alone.
        sim = make([host(1., .5, 0., 1., omega=2.)], {0: RING})
        sim.step_collision()
        self.assertEqual(sim.events[-1].kind, "specular")
        np.testing.assert_allclose(sim.world.bodies.vel[0], [0., -1.])
        self.assertEqual(sim.world.bodies.omega[0], 2.)
        # Post: rough, energy kept, momentum change equals the support impulse.
        sim = make([host(.5, .55, 1., 0., omega=1.)], {0: RING},
                   posts=[CircularPost((1., .5), .05)])
        e0, p0, _ = invariants(sim)
        sim.step_collision()
        self.assertEqual(sim.events[-1].kind, "post")
        e1, p1, _ = invariants(sim)
        self.assertAlmostEqual(e1, e0, delta=1e-13)
        np.testing.assert_allclose(p1 - p0, -sim.post_impulses()[0], atol=1e-13)


@unittest.skipUnless(edmd.numba_available(), "Numba optional dependency not installed")
class RingDynamicsTests(unittest.TestCase):
    def test_closed_box_with_hosts_conserves_energy_and_never_overlaps(self):
        rng = np.random.default_rng(4)
        specs, rings = gas_with_hosts(rng, 6, 200, 3., 1.5)
        sim = make(specs, rings, 3., 1.5)
        e0 = sim.energy()
        for k in range(1, 41):
            sim.advance_to(k*.5)
        kinds = {}
        for e in sim.events:
            key = e.metadata["feature"] if e.kind == "host" else e.kind
            kinds[key] = kinds.get(key, 0) + 1
        for feature in ("outer", "inner", "host_host"):
            self.assertGreater(kinds.get(feature, 0), 0, kinds)
        self.assertGreater(kinds.get("cap_plus", 0) + kinds.get("cap_minus", 0), 0, kinds)
        self.assertLessEqual(abs(sim.energy() - e0), 1e-11*e0)
        self.assertLessEqual(sim.max_penetration, 1e-10)

    def test_host_registration_agrees_with_a_single_cell(self):
        rng = np.random.default_rng(8)
        specs, rings = gas_with_hosts(rng, 4, 80, 2., 1.)
        grid = make(specs, rings)
        one = make(specs, rings)
        one.fp[edmd.FP_CX] *= one.ip[0]
        one.fp[edmd.FP_CY] *= one.ip[1]
        one.ip[:] = 1
        one.head = np.full(1, -1, np.int64)
        one.cps, one.cpi = np.zeros(2, np.int64), np.zeros(0, np.int64)
        one.reg = np.zeros((1, 4), np.int8)
        edmd.build_grid(one._bodies(), one._grid(), one._hosts(), 0.)
        edmd.register_all(one._bodies(), one._grid(), one._hosts(), 0.)
        one._rebuild()
        grid.advance_to(2.)
        one.advance_to(2.)
        pairs = list(zip(grid.events, one.events))[:200]
        self.assertEqual(len(pairs), 200)
        self.assertGreater(sum(a.kind == "host" for a, _ in pairs), 10)
        for a, b in pairs:
            self.assertEqual((a.kind, a.participants, a.metadata),
                             (b.kind, b.participants, b.metadata))
            self.assertAlmostEqual(b.time, a.time, delta=1e-8)

    def test_partial_roughness_is_reproducible_through_checkpoints(self):
        rng = np.random.default_rng(9)
        specs, rings = gas_with_hosts(rng, 4, 100, 2., 1.)
        sim = make(specs, rings, thermal=BoundaryKind.HOT, roughness=.5)
        sim.advance_to(1.)
        cp = sim.checkpoint()
        first = sim.advance_to(3.)
        events = [(e.time, e.kind, e.participants) for e in sim.events]
        rough = [e.metadata["rough"] for e in sim.events if e.kind == "host"]
        self.assertTrue(any(rough) and not all(rough))
        sim.restore(cp)
        second = sim.advance_to(3.)
        np.testing.assert_array_equal(first.position, second.position)
        np.testing.assert_array_equal(first.angle, second.angle)
        self.assertEqual([(e.time, e.kind, e.participants) for e in sim.events], events)

    def test_discs_enter_and_leave_the_cavities(self):
        rng = np.random.default_rng(12)
        specs, rings = gas_with_hosts(rng, 4, 150, 2., 1.)
        sim = make(specs, rings, thermal=BoundaryKind.HOT)
        inside_counts = []
        for k in range(1, 201):
            snap = sim.advance_to(k*.5)
            p = snap.position
            inside = 0
            for h in range(4):
                inside += int(np.sum(np.linalg.norm(p[4:] - p[h], axis=1) < RING.inner_radius))
            inside_counts.append(inside)
        self.assertGreater(max(inside_counts), 0)
        self.assertGreater(len(set(inside_counts)), 2)     # it fluctuates
        self.assertLessEqual(abs(sim.snapshot().energy_residual), 1e-10*sim.energy())

    def test_equipartition_between_discs_hosts_spin_and_rotation(self):
        for roughness in (1., .3):
            rng = np.random.default_rng(31)
            specs, rings = gas_with_hosts(rng, 6, 200, 3., 1.5, temperature=.4)
            world = World(BodyArrays.from_specs(specs),
                          [w for w in box_walls(3., 1.5, BoundaryKind.HOT, 1.)],
                          rings=rings, contact_roughness=roughness)
            for w in world.walls:      # every wall thermal at T = 1
                w.kind, w.temperature = BoundaryKind.HOT, 1.
            sim = EdmdSimulation(world, record_events=False)
            sim.advance_to(200.)                         # relax from T = 0.4
            s = sim.world.bodies
            hosts, discs = slice(0, 6), slice(6, None)
            sums = np.zeros(4)
            samples = 0
            for k in range(1, 1001):
                sim.advance_to(200. + k)
                sums += [np.mean(s.mass[discs]*np.sum(s.vel[discs]**2, 1))/2,
                         np.mean(s.inertia[discs]*s.omega[discs]**2),
                         np.mean(s.mass[hosts]*np.sum(s.vel[hosts]**2, 1))/2,
                         np.mean(s.inertia[hosts]*s.omega[hosts]**2)]
                samples += 1
            disc_t, disc_spin, host_t, host_rot = sums/samples
            with self.subTest(roughness=roughness):
                self.assertAlmostEqual(disc_t, 1., delta=.03)
                self.assertAlmostEqual(disc_spin, 1., delta=.06)
                self.assertAlmostEqual(host_t, 1., delta=.12)
                self.assertAlmostEqual(host_rot, 1., delta=.15)


if __name__ == "__main__":
    unittest.main()
