from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math


@dataclass(frozen=True)
class RunConfig:
    preset: str = "gas_box"
    seed: int = 123
    duration: float = 10.0
    particles: int = 32
    temperature: float = 1.0
    sample_interval: float = 0.05
    max_horizon: float = 0.05
    cycles: int | None = None
    reversed_cycle: bool = False
    shaft_mode: str = "controlled"
    shaft_speed: float = 0.15
    transient_cycles: int = 2
    efficiency_min_cycles: int = 8
    pair_search: str = "grid"
    numeric_backend: str = "auto"
    wall_search: str = "bounded"
    wall_kernel: str = "auto"
    penetration_kernel: str = "auto"

    def validate(self) -> None:
        if self.duration <= 0 or self.particles <= 0 or self.temperature <= 0:
            raise ValueError("duration, particles, and temperature must be positive")
        if self.sample_interval <= 0 or self.max_horizon <= 0:
            raise ValueError("sampling and horizon intervals must be positive")
        if self.shaft_mode not in ("controlled", "free"):
            raise ValueError("shaft_mode must be 'controlled' or 'free'")
        if not math.isfinite(self.shaft_speed) or self.shaft_speed <= 0:
            raise ValueError("shaft_speed must be finite and positive")
        if self.transient_cycles < 0 or self.efficiency_min_cycles < 4:
            raise ValueError("require nonnegative transients and at least four efficiency cycles")
        if self.pair_search not in ("grid", "all"):
            raise ValueError("pair_search must be grid or all")
        if self.numeric_backend not in ("auto", "python", "numba"):
            raise ValueError("numeric_backend must be auto, python, or numba")
        if self.wall_search not in ("bounded", "all"):
            raise ValueError("wall_search must be bounded or all")
        if self.wall_kernel not in ("auto", "python"):
            raise ValueError("wall_kernel must be auto or python")
        if self.penetration_kernel not in ("auto", "python"):
            raise ValueError("penetration_kernel must be auto or python")

    def digest(self) -> str:
        return hashlib.sha256(json.dumps(asdict(self), sort_keys=True).encode()).hexdigest()
