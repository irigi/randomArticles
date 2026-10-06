from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math
import numpy as np

from .contacts import apply_impulse, inverse_effective_mass, relative_normal_velocity
from .events import Contact
from .state import BodyArrays


class BoundaryKind(str, Enum):
    SPECULAR = "specular"
    HOT = "hot"
    COLD = "cold"


@dataclass
class Wall:
    point: np.ndarray
    inward_normal: np.ndarray
    length: float
    kind: BoundaryKind = BoundaryKind.SPECULAR
    temperature: float | None = None
    name: str = "wall"
    velocity: np.ndarray | None = None

    def __post_init__(self):
        self.point = np.asarray(self.point, dtype=float)
        self.inward_normal = np.asarray(self.inward_normal, dtype=float)
        self.inward_normal /= np.linalg.norm(self.inward_normal)
        self.velocity = np.zeros(2) if self.velocity is None else np.asarray(self.velocity, float)
        if self.length <= 0:
            raise ValueError("wall length must be positive")
        if self.kind != BoundaryKind.SPECULAR and (self.temperature is None or self.temperature <= 0):
            raise ValueError("thermal wall requires positive temperature")

    def point_at(self, time: float) -> np.ndarray:
        return self.point + self.velocity*time

    def velocity_at(self, time: float) -> np.ndarray:
        return self.velocity

    def speed_bound(self, t0: float, t1: float) -> float:
        return float(np.linalg.norm(self.velocity))

    def kind_at(self, time: float) -> BoundaryKind:
        return self.kind

    def temperature_at(self, time: float) -> float | None:
        return self.temperature


@dataclass
class CircularPost:
    """Fixed, specular circular obstacle, e.g. one post of a porous membrane."""

    center: np.ndarray
    radius: float
    name: str = "post"

    def __post_init__(self):
        self.center = np.asarray(self.center, dtype=float).reshape(2)
        if not (math.isfinite(self.radius) and self.radius > 0):
            raise ValueError("post radius must be positive")


class SegmentWall(Wall):
    """Stationary finite two-sided wall with explicit material thickness."""

    def __init__(self, start, end, thickness: float = 0.01,
                 kind: BoundaryKind = BoundaryKind.SPECULAR,
                 temperature: float | None = None, name: str = "segment"):
        self.start = np.asarray(start, float); self.end = np.asarray(end, float)
        edge = self.end-self.start
        length = float(np.linalg.norm(edge))
        if length <= 0 or thickness < 0:
            raise ValueError("segment length must be positive and thickness nonnegative")
        self.thickness = thickness
        normal = np.array([-edge[1], edge[0]])/length
        super().__init__(.5*(self.start+self.end), normal, length, kind, temperature, name)


def resolve_wall(state: BodyArrays, body: int, wall: Wall, point: np.ndarray,
                 rng: np.random.Generator, time: float = 0.0,
                 contact_normal: np.ndarray | None = None,
                 diagnostics: dict | None = None) -> tuple[float, float, float]:
    """Return scalar impulse, heat into body, and work on body."""
    # Contact normal points from body toward the wall.
    n = -wall.inward_normal if contact_normal is None else np.asarray(contact_normal, float)
    c = Contact(body, None, tuple(point), tuple(n))
    before = state.kinetic_energy()
    wall_velocity = wall.velocity_at(time)
    kind = wall.kind_at(time)
    temperature = wall.temperature_at(time)
    g = relative_normal_velocity(state, c) + float(wall_velocity @ n)
    d = inverse_effective_mass(state, c)
    if g >= 0 or d <= 0:
        if diagnostics is not None and kind != BoundaryKind.SPECULAR:
            mode_energy=float(g*g/(2*d)) if d > 0 else 0.0
            diagnostics.update(
                reservoir_temperature=float(temperature),
                incoming_normal_mode_energy=mode_energy,
                outgoing_normal_mode_energy=mode_energy)
        return 0.0, 0.0, 0.0
    if kind == BoundaryKind.SPECULAR:
        target = -g
    else:
        # Flux-weighted contact-normal equilibrium mode.
        target = math.sqrt(-2.0 * temperature * d * math.log(max(rng.random(), np.finfo(float).tiny)))
    impulse = (target-g)/d
    apply_impulse(state, c, impulse)
    after = state.kinetic_energy()
    delta = after-before
    # The impulse on the body is -j*n, so prescribed-boundary work on the
    # body is that impulse dotted with the boundary velocity.
    work = float(-impulse * (wall_velocity @ n))
    heat = delta-work if kind != BoundaryKind.SPECULAR else 0.0
    if diagnostics is not None and kind != BoundaryKind.SPECULAR:
        diagnostics.update(
            reservoir_temperature=float(temperature),
            incoming_normal_mode_energy=float(g*g/(2*d)),
            outgoing_normal_mode_energy=float(target*target/(2*d)))
    return impulse, heat, work
