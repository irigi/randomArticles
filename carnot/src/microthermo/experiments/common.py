from __future__ import annotations

import math
import numpy as np
from ..core.boundaries import BoundaryKind, Wall
from ..core.state import BodySpec, Shape


def box_walls(width: float, height: float, thermal: BoundaryKind | None = None,
              temperature: float = 1.0) -> list[Wall]:
    k = thermal or BoundaryKind.SPECULAR
    temp = temperature if thermal else None
    return [
        Wall(np.array([0., height/2]), np.array([1., 0.]), height, k, temp, "left"),
        Wall(np.array([width, height/2]), np.array([-1., 0.]), height, BoundaryKind.SPECULAR, None, "right"),
        Wall(np.array([width/2, 0.]), np.array([0., 1.]), width, BoundaryKind.SPECULAR, None, "bottom"),
        Wall(np.array([width/2, height]), np.array([0., -1.]), width, BoundaryKind.SPECULAR, None, "top"),
    ]


def particle_specs(count: int, width: float, height: float, radius: float,
                   temperature: float, rng: np.random.Generator,
                   shape: Shape = Shape.DISC, x_max: float | None = None) -> list[BodySpec]:
    xmax = width if x_max is None else x_max
    specs: list[BodySpec] = []
    max_tries = 10000*count
    tries = 0
    while len(specs) < count and tries < max_tries:
        tries += 1
        p = np.array([rng.uniform(radius, xmax-radius), rng.uniform(radius, height-radius)])
        if any(np.linalg.norm(p-np.asarray(s.position)) < 2.05*radius for s in specs):
            continue
        v = rng.normal(0., math.sqrt(temperature), 2)
        omega = rng.normal(0., math.sqrt(temperature/(0.25*radius*radius))) if shape != Shape.DISC else 0.0
        specs.append(BodySpec(tuple(p), tuple(v), radius=radius, shape=shape, omega=omega))
    if len(specs) != count:
        raise ValueError("could not place particles without overlap; lower count or radius")
    # Remove center-of-mass drift without changing the requested temperature materially.
    mean = np.mean([s.velocity for s in specs], axis=0)
    for s in specs:
        s.velocity = tuple(np.asarray(s.velocity)-mean)
    return specs

