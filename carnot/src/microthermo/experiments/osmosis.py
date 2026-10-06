"""Osmosis lab: ring hosts with binding wells behind a disc-permeable membrane.

The box is split by a column of posts at its middle. Small discs pass the
gaps; ring hosts are larger than the gaps and stay on the left. The right
chamber holds only discs, so its density gives the discs' chemical
potential directly; on the left, discs are free or bound in the hosts'
wells. All four outer walls are thermal at the run temperature.

Sizing (see docs/PLAN_OSMOSIS.md §2): the hosts' outer circles cover a fixed
fraction of the left chamber, so they shrink as their number grows, while
wall thickness, mouth width and well clearance are tied to the disc size.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from ..config import RunConfig
from ..core.boundaries import BoundaryKind
from ..core.state import BodyArrays, BodySpec, RingGeometry, Shape
from ..runner.simulation import World
from .common import box_walls
from .membrane import Membrane, build_membrane

WIDTH, HEIGHT = 4.0, 2.0
DISC_RADIUS = 0.02
WALL_THICKNESS = 0.03
POST_RADIUS = 0.07
POST_GAP = 0.10


@dataclass(frozen=True)
class OsmosisDesign:
    """Geometry and the ideal-point reference for one osmosis configuration."""

    hosts: int
    ring: RingGeometry
    membrane: Membrane
    disc_radius: float
    host_mass: float
    temperature: float
    # Areas accessible to disc centres (Monte Carlo, empty box with membrane).
    right_area: float
    left_area: float
    # Area accessible to host centres in the left chamber.
    host_area: float

    @property
    def binding_energy(self) -> float:
        return self.ring.well_depth

    @property
    def well_area(self) -> float:
        return math.pi*self.ring.well_radius**2

    @property
    def free_left_area(self) -> float:
        """Left area open to free disc centres: hosts' footprints removed,
        their cavities outside the wells added back."""
        r, g = self.disc_radius, self.ring
        per_host = (math.pi*(g.outer_radius + r)**2 - math.pi*(g.inner_radius - r)**2
                    + self.well_area)
        return self.left_area - self.hosts*per_host

    def reference(self, discs: int, temperature: float | None = None) -> dict[str, float]:
        """Ideal non-interacting discs, total number fixed (canonical)."""
        t = self.temperature if temperature is None else temperature
        weight = self.well_area*math.exp(self.binding_energy/t)
        c = discs/(self.right_area + self.free_left_area + self.hosts*weight)
        density = self.hosts/self.host_area
        # Hard-disc second virial coefficient for the hosts' outer circles.
        b2 = math.pi*(2*self.ring.outer_radius)**2/2
        return {"concentration": c, "right": c*self.right_area,
                "left": discs - c*self.right_area,
                "free_left": c*self.free_left_area,
                "bound": c*self.hosts*weight, "per_host": c*weight,
                "osmotic_pressure_ideal": density*t,
                "osmotic_pressure_virial": density*t*(1 + b2*density)}


def _accessible_areas(membrane: Membrane, radius: float, samples: int = 400000,
                      seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    p = rng.uniform((radius, radius), (WIDTH - radius, HEIGHT - radius), (samples, 2))
    clear = np.ones(samples, bool)
    for post in membrane.posts:
        clear &= np.hypot(*(p - post.center).T) > post.radius + radius
    box = (WIDTH - 2*radius)*(HEIGHT - 2*radius)
    left = clear & (p[:, 0] < membrane.x)
    right = clear & (p[:, 0] >= membrane.x)
    return box*np.mean(right), box*np.mean(left)


def design_osmosis(config: RunConfig) -> OsmosisDesign:
    hard = config.preset == "osmosis_hard"
    depth = 0.0 if hard else (1.5 if config.binding_energy is None else config.binding_energy)
    r = DISC_RADIUS
    mouth = 6*r if config.mouth_width is None else config.mouth_width
    membrane_x = WIDTH/2
    left_area = membrane_x*HEIGHT
    outer = math.sqrt(config.host_area_fraction*left_area/(math.pi*config.hosts))
    inner = outer - WALL_THICKNESS
    well = inner - 2*r
    if well < 4*r:
        raise ValueError(
            f"{config.hosts} hosts at area fraction {config.host_area_fraction:g} leave a well "
            f"of radius {well:.3g}, below two disc diameters ({4*r:g}); use fewer hosts "
            "or a larger area fraction")
    ring = RingGeometry(inner, outer, mouth, well_radius=well, well_depth=depth)
    membrane = build_membrane(membrane_x, 0., HEIGHT, POST_RADIUS, POST_GAP, r,
                              blocked_radius=outer, mouth_width=mouth)
    right, left = _accessible_areas(membrane, r)
    host_area = ((membrane_x - POST_RADIUS - 2*outer)*(HEIGHT - 2*outer))
    return OsmosisDesign(config.hosts, ring, membrane, r, config.host_mass,
                         config.temperature, right, left, host_area)


class OsmosisExperiment:
    def build(self, config: RunConfig, rng: np.random.Generator) -> World:
        design = design_osmosis(config)
        ring, membrane, r = design.ring, design.membrane, design.disc_radius
        t, m_host = config.temperature, config.host_mass
        clearance = .01
        host_inertia = ring.uniform_inertia(m_host)
        x_max = membrane.x - membrane.post_radius - ring.outer_radius - clearance
        centres: list[np.ndarray] = []
        for _ in range(100000):
            if len(centres) == design.hosts:
                break
            p = rng.uniform((ring.outer_radius + clearance,)*2,
                            (x_max, HEIGHT - ring.outer_radius - clearance))
            if all(np.hypot(*(p - q)) > 2*ring.outer_radius + clearance for q in centres):
                centres.append(p)
        else:
            raise ValueError("could not place the hosts without overlap")
        specs = [BodySpec(tuple(p), tuple(rng.normal(0, math.sqrt(t/m_host), 2)),
                          mass=m_host, radius=ring.outer_radius, shape=Shape.HOST,
                          angle=float(rng.uniform(-math.pi, math.pi)),
                          omega=float(rng.normal(0, math.sqrt(t/host_inertia))),
                          inertia=host_inertia)
                 for p in centres]
        disc_inertia = .5*r*r
        x_low = (membrane.x + membrane.post_radius + r + clearance
                 if config.discs_start == "right" else r + clearance)
        placed: list[np.ndarray] = []
        for _ in range(1000*config.particles):
            if len(placed) == config.particles:
                break
            p = rng.uniform((x_low, r + clearance), (WIDTH - r - clearance, HEIGHT - r - clearance))
            if any(np.hypot(*(p - q)) < ring.outer_radius + r + clearance for q in centres):
                continue
            if any(np.hypot(*(p - post.center)) < post.radius + r + clearance
                   for post in membrane.posts):
                continue
            if placed and np.min(np.hypot(*(np.asarray(placed) - p).T)) < 2*r + clearance:
                continue
            placed.append(p)
        else:
            raise ValueError("could not place the discs without overlap; lower the count")
        for p in placed:
            specs.append(BodySpec(tuple(p), tuple(rng.normal(0, math.sqrt(t), 2)), radius=r,
                                  omega=float(rng.normal(0, math.sqrt(t/disc_inertia))),
                                  inertia=disc_inertia))
        walls = box_walls(WIDTH, HEIGHT)
        for wall in walls:
            wall.kind, wall.temperature = BoundaryKind.HOT, t
        reference = design.reference(config.particles)
        metadata = {
            "name": config.preset, "membrane_x": membrane.x, "membrane_gap": membrane.gap,
            "post_radius": membrane.post_radius, "hosts": design.hosts,
            "ring": {"inner_radius": ring.inner_radius, "outer_radius": ring.outer_radius,
                     "mouth_width": ring.mouth_width, "well_radius": ring.well_radius,
                     "well_depth": ring.well_depth},
            "disc_radius": r, "host_mass": m_host, "temperature": t,
            "contact_roughness": config.rough_fraction,
            "accessible_area_right": design.right_area,
            "accessible_area_left": design.left_area,
            "free_left_area": design.free_left_area, "well_area": design.well_area,
            "host_accessible_area": design.host_area, "reference": reference,
            "interpretation": "geometric and energetic association behind a "
                              "disc-permeable membrane; hosts meet walls, posts and "
                              "each other as their outer circles",
        }
        return World(BodyArrays.from_specs(specs), walls, posts=list(membrane.posts),
                     rings={i: ring for i in range(design.hosts)},
                     contact_roughness=config.rough_fraction, metadata=metadata)
