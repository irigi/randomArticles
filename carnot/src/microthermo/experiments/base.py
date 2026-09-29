from __future__ import annotations

from abc import ABC, abstractmethod
import numpy as np
from ..config import RunConfig
from ..runner.simulation import World


class Experiment(ABC):
    @abstractmethod
    def build(self, config: RunConfig, rng: np.random.Generator) -> World: ...

    def instruments(self) -> list[str]:
        return ["energy", "temperature", "pressure"]

    def guide(self) -> list[str]:
        return []

