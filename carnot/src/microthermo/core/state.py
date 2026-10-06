from __future__ import annotations

from dataclasses import dataclass, field
from enum import IntEnum
import numpy as np


class Shape(IntEnum):
    DISC = 0
    TRIANGLE = 1
    HOST = 2


@dataclass
class BodySpec:
    position: tuple[float, float]
    velocity: tuple[float, float]
    mass: float = 1.0
    radius: float = 0.05
    shape: Shape = Shape.DISC
    angle: float = 0.0
    omega: float = 0.0
    inertia: float | None = None
    species: int = 0
    dynamic: bool = True
    vertices: list[tuple[float, float]] | None = None


@dataclass
class BodyArrays:
    """Structure-of-arrays state used by the event loop."""

    pos: np.ndarray
    vel: np.ndarray
    angle: np.ndarray
    omega: np.ndarray
    mass: np.ndarray
    inertia: np.ndarray
    radius: np.ndarray
    shape: np.ndarray
    species: np.ndarray
    dynamic: np.ndarray
    version: np.ndarray
    ids: np.ndarray
    polygons: dict[int, np.ndarray] = field(default_factory=dict)
    packed_triangles: np.ndarray | None = None

    @classmethod
    def from_specs(cls, specs: list[BodySpec]) -> "BodyArrays":
        if not specs:
            raise ValueError("at least one body is required")
        n = len(specs)
        pos = np.asarray([s.position for s in specs], dtype=np.float64)
        vel = np.asarray([s.velocity for s in specs], dtype=np.float64)
        mass = np.asarray([s.mass for s in specs], dtype=np.float64)
        radius = np.asarray([s.radius for s in specs], dtype=np.float64)
        shape = np.asarray([int(s.shape) for s in specs], dtype=np.int16)
        dynamic = np.asarray([s.dynamic for s in specs], dtype=bool)
        if np.any(mass <= 0) or np.any(radius <= 0):
            raise ValueError("mass and circumradius must be positive")
        inertia = np.empty(n, dtype=np.float64)
        polygons: dict[int, np.ndarray] = {}
        for i, s in enumerate(specs):
            if s.inertia is not None:
                inertia[i] = s.inertia
            elif s.shape == Shape.DISC:
                inertia[i] = np.inf  # smooth discs have no spin degree of freedom
            elif s.shape == Shape.TRIANGLE:
                side = np.sqrt(3.0) * s.radius
                inertia[i] = s.mass * side * side / 12.0
            else:
                inertia[i] = 0.5 * s.mass * s.radius * s.radius
            if s.vertices is not None:
                v = np.asarray(s.vertices, dtype=np.float64)
                if v.ndim != 2 or v.shape[1] != 2 or len(v) < 3:
                    raise ValueError("polygon vertices must be an Nx2 array")
                polygons[i] = v
            elif s.shape == Shape.TRIANGLE:
                a = np.array([np.pi / 2, np.pi / 2 + 2*np.pi/3, np.pi / 2 + 4*np.pi/3])
                polygons[i] = s.radius * np.column_stack((np.cos(a), np.sin(a)))
        if np.any((inertia <= 0) & np.isfinite(inertia)):
            raise ValueError("finite inertia must be positive")
        packed_triangles = None
        if len(polygons) == n and all(len(polygons[i]) == 3 for i in range(n)):
            packed_triangles = np.stack([polygons[i] for i in range(n)])
            polygons = {i: packed_triangles[i] for i in range(n)}
        return cls(pos, vel, np.asarray([s.angle for s in specs], float),
                   np.asarray([s.omega for s in specs], float), mass, inertia,
                   radius, shape, np.asarray([s.species for s in specs], np.int16),
                   dynamic, np.zeros(n, np.int64), np.arange(n, dtype=np.int64),
                   polygons, packed_triangles)

    def copy(self) -> "BodyArrays":
        if self.packed_triangles is not None:
            packed = self.packed_triangles.copy()
            polygons = {i: packed[i] for i in range(self.n)}
        else:
            packed = None
            polygons = {k: v.copy() for k, v in self.polygons.items()}
        return BodyArrays(*(getattr(self, name).copy() for name in
            ("pos", "vel", "angle", "omega", "mass", "inertia", "radius",
             "shape", "species", "dynamic", "version", "ids")),
            polygons, packed)

    @property
    def n(self) -> int:
        return len(self.mass)

    def kinetic_energy(self) -> float:
        linear = 0.5 * np.sum(self.mass * np.sum(self.vel * self.vel, axis=1))
        rotational = 0.5 * np.sum(np.where(np.isfinite(self.inertia), self.inertia, 0.0) * self.omega**2)
        return float(linear + rotational)

    def momentum(self) -> np.ndarray:
        return np.sum(self.mass[:, None] * self.vel, axis=0)

    def angular_momentum(self) -> float:
        orbital = np.sum(self.pos[:, 0] * self.mass * self.vel[:, 1] - self.pos[:, 1] * self.mass * self.vel[:, 0])
        spin = np.sum(np.where(np.isfinite(self.inertia), self.inertia, 0.0) * self.omega)
        return float(orbital + spin)

    def advance(self, dt: float) -> None:
        moving = self.dynamic
        self.pos[moving] += self.vel[moving] * dt
        self.angle[moving] = np.remainder(self.angle[moving] + self.omega[moving] * dt + np.pi, 2*np.pi) - np.pi


@dataclass(frozen=True)
class RingGeometry:
    """Ring host: an annular wall with one mouth and rounded mouth ends.

    The solid is every point within ``cap`` of an arc of radius ``mid``;
    the arc leaves out the mouth, centred on body angle 0. ``mouth_width``
    is the clear chord between the two rounded ends.
    """

    inner_radius: float
    outer_radius: float
    mouth_width: float

    def __post_init__(self) -> None:
        if not (0 < self.inner_radius < self.outer_radius and self.mouth_width > 0):
            raise ValueError("ring needs 0 < inner radius < outer radius and a positive mouth")
        if self.mouth_width + 2*self.cap >= 2*self.mid:
            raise ValueError("ring mouth is wider than the ring")

    @property
    def mid(self) -> float:
        return 0.5*(self.inner_radius + self.outer_radius)

    @property
    def cap(self) -> float:
        return 0.5*(self.outer_radius - self.inner_radius)

    @property
    def mouth_half_angle(self) -> float:
        """Angle of each rounded end's centre from the mouth axis."""
        return float(np.arcsin((self.mouth_width + 2*self.cap)/(2*self.mid)))

    def uniform_inertia(self, mass: float) -> float:
        """Moment of inertia of a full uniform annulus (the declared balanced mass)."""
        return 0.5*mass*(self.inner_radius**2 + self.outer_radius**2)


@dataclass(frozen=True)
class Tolerances:
    geometry: float = 1e-10
    time: float = 1e-12
    velocity: float = 1e-12
    energy: float = 1e-10
