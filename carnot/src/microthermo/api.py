from __future__ import annotations

import numpy as np
from .config import RunConfig
from .experiments import build_preset
from .runner.simulation import Simulation, Snapshot


def load_preset(config: RunConfig) -> Simulation:
    config.validate()
    world = build_preset(config, np.random.default_rng(config.seed))
    return Simulation(world, config.seed, config.max_horizon,
                      pair_search=config.pair_search,
                      numeric_backend=config.numeric_backend,
                      wall_search=config.wall_search,
                      wall_kernel=config.wall_kernel,
                      penetration_kernel=config.penetration_kernel)


__all__ = ["Simulation", "Snapshot", "RunConfig", "load_preset"]
