from __future__ import annotations

import numpy as np
from .common import box_walls, particle_specs
from .carnot import CarnotExperiment
from .labyrinth import LabyrinthExperiment
from .permeability import PermeabilityExperiment
from ..config import RunConfig
from ..core.boundaries import BoundaryKind
from ..core.state import BodyArrays, Shape
from ..runner.simulation import World


def preset_names() -> tuple[str, ...]:
    return ("gas_box", "triangle_equipartition", "carnot_discs", "carnot_triangles",
            "labyrinth", "labyrinth_energetic", "selective_membrane")


def build_preset(config: RunConfig, rng: np.random.Generator) -> World:
    if config.preset == "carnot_discs": return CarnotExperiment(False).build(config, rng)
    if config.preset == "carnot_triangles": return CarnotExperiment(True).build(config, rng)
    if config.preset == "labyrinth": return LabyrinthExperiment().build(config, rng)
    if config.preset == "labyrinth_energetic": return LabyrinthExperiment(config.temperature).build(config, rng)
    if config.preset == "selective_membrane": return PermeabilityExperiment().build(config, rng)
    if config.preset in ("gas_box", "triangle_equipartition"):
        shape = Shape.TRIANGLE if config.preset == "triangle_equipartition" else Shape.DISC
        specs = particle_specs(config.particles, 2., 1., 0.025, config.temperature, rng, shape)
        if shape == Shape.TRIANGLE:
            for s in specs:
                s.omega *= 0.1
        return World(BodyArrays.from_specs(specs), box_walls(2., 1., BoundaryKind.HOT,
                     config.temperature), metadata={"name": config.preset})
    raise ValueError(f"unknown preset {config.preset!r}; choose from {', '.join(preset_names())}")
