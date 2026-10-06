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
    particles: int | None = None
    temperature: float = 1.0
    sample_interval: float = 0.05
    max_horizon: float = 0.05
    cycles: int | None = None
    reversed_cycle: bool = False
    shaft_mode: str = "controlled"
    shaft_speed: float = 0.15
    cold_jacket: bool = False
    hot_jacket: bool = False
    cam_fractions: tuple[float, float, float, float] = (.25, .25, .25, .25)
    transient_cycles: int = 2
    efficiency_min_cycles: int = 8
    pair_search: str = "grid"
    numeric_backend: str = "auto"
    wall_search: str = "bounded"
    wall_kernel: str = "auto"
    penetration_kernel: str = "auto"
    pair_kernel: str = "auto"
    cam_kernel: str = "auto"
    carnot_radius_scale: float = 1.0
    # Carnot only: gas starts at this temperature instead of `temperature`,
    # which still sets the reservoirs (T_hot = 1.5*temperature).
    initial_temperature: float | None = None
    # Carnot only: T_hot/T_cold. The cam adiabats are recalibrated to it.
    temperature_ratio: float = 2.0
    # "reference" (bounded-horizon scheduler) or "edmd" (compiled
    # event-driven kernel for smooth discs in a wall box).
    engine: str = "reference"

    def __post_init__(self) -> None:
        object.__setattr__(self, "cam_fractions", tuple(self.cam_fractions))
        if self.particles is None:
            object.__setattr__(self, "particles",
                               96 if self.preset == "carnot_triangles" else 32)

    def validate(self) -> None:
        if self.duration <= 0 or self.particles <= 0 or self.temperature <= 0:
            raise ValueError("duration, particles, and temperature must be positive")
        if self.sample_interval <= 0 or self.max_horizon <= 0:
            raise ValueError("sampling and horizon intervals must be positive")
        if self.shaft_mode not in ("controlled", "free"):
            raise ValueError("shaft_mode must be 'controlled' or 'free'")
        if not math.isfinite(self.shaft_speed) or self.shaft_speed <= 0:
            raise ValueError("shaft_speed must be finite and positive")
        if self.cold_jacket and not self.preset.startswith("carnot_"):
            raise ValueError("cold jacket applies only to Carnot presets")
        if self.hot_jacket and not self.preset.startswith("carnot_"):
            raise ValueError("hot jacket applies only to Carnot presets")
        if not math.isfinite(self.carnot_radius_scale) or self.carnot_radius_scale <= 0:
            raise ValueError("carnot_radius_scale must be finite and positive")
        if self.carnot_radius_scale != 1.0 and not self.preset.startswith("carnot_"):
            raise ValueError("carnot_radius_scale applies only to Carnot presets")
        if len(self.cam_fractions) != 4 or any(
                not math.isfinite(x) or x <= 0 for x in self.cam_fractions):
            raise ValueError("cam_fractions must contain four positive finite values")
        if not math.isclose(sum(self.cam_fractions), 1., rel_tol=0., abs_tol=1e-12):
            raise ValueError("cam_fractions must sum to one")
        if self.cam_fractions != (.25, .25, .25, .25) and not self.preset.startswith("carnot_"):
            raise ValueError("cam_fractions applies only to Carnot presets")
        if self.initial_temperature is not None and not (
                math.isfinite(self.initial_temperature) and self.initial_temperature > 0):
            raise ValueError("initial_temperature must be finite and positive")
        if not math.isfinite(self.temperature_ratio) or self.temperature_ratio <= 1:
            raise ValueError("temperature_ratio must be finite and above one")
        if ((self.initial_temperature is not None or self.temperature_ratio != 2.0)
                and not self.preset.startswith("carnot_")):
            raise ValueError("initial_temperature and temperature_ratio apply only to Carnot presets")
        if self.transient_cycles < 0 or self.efficiency_min_cycles < 4:
            raise ValueError("require nonnegative transients and at least four efficiency cycles")
        if self.pair_search not in ("grid", "sweep", "all"):
            raise ValueError("pair_search must be grid, sweep, or all")
        if self.numeric_backend not in ("auto", "python", "numba"):
            raise ValueError("numeric_backend must be auto, python, or numba")
        if self.wall_search not in ("bounded", "all"):
            raise ValueError("wall_search must be bounded or all")
        if self.wall_kernel not in ("auto", "python"):
            raise ValueError("wall_kernel must be auto or python")
        if self.penetration_kernel not in ("auto", "python"):
            raise ValueError("penetration_kernel must be auto or python")
        if self.pair_kernel not in ("auto", "scalar"):
            raise ValueError("pair_kernel must be auto or scalar")
        if self.cam_kernel not in ("auto", "python"):
            raise ValueError("cam_kernel must be auto or python")
        if self.engine not in ("reference", "edmd"):
            raise ValueError("engine must be reference or edmd")

    def digest(self) -> str:
        data = asdict(self)
        # Fields added later are omitted at their defaults so older run
        # hashes, archives, and checkpoints remain valid.
        for key, default in (("initial_temperature", None), ("temperature_ratio", 2.0),
                             ("engine", "reference")):
            if data[key] == default:
                del data[key]
        return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()
