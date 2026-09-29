from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ApparatusComponent:
    """Immutable render and instrument state for one named mechanism part.

    Points use the gas cylinder's physical coordinate system. A closed path
    repeats its first point. Energy belongs to this component only; zero-energy
    cam and selector geometry does not add another mechanical degree of freedom.
    """

    name: str
    kind: str
    points: tuple[tuple[float, float], ...]
    kinetic_energy: float = 0.0
    potential_energy: float = 0.0
    work_output: float = 0.0
    velocity: tuple[float, float] = (0.0, 0.0)
    angular_velocity: float = 0.0
    state: str = ""
