from __future__ import annotations

import numpy as np
from .config import RunConfig
from .experiments import build_preset
from .runner.edmd_simulation import EdmdSimulation
from .runner.simulation import Simulation, Snapshot


def load_preset(config: RunConfig) -> Simulation | EdmdSimulation:
    config.validate()
    world = build_preset(config, np.random.default_rng(config.seed))
    if config.engine == "edmd":
        return EdmdSimulation(world, config.seed)
    return Simulation(world, config.seed, config.max_horizon,
                      pair_search=config.pair_search,
                      numeric_backend=config.numeric_backend,
                      wall_search=config.wall_search,
                      wall_kernel=config.wall_kernel,
                      penetration_kernel=config.penetration_kernel,
                      pair_kernel=config.pair_kernel,
                      cam_kernel=config.cam_kernel)


__all__ = ["EdmdSimulation", "Simulation", "Snapshot", "RunConfig", "load_preset"]
