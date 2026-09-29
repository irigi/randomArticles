"""Analytical point-gas Carnot comparison for the designed cam geometry."""

from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np


@dataclass(frozen=True)
class IdealBranch:
    name: str
    area: tuple[float, ...]
    pressure: tuple[float, ...]


@dataclass(frozen=True)
class IdealCarnotReference:
    branches: tuple[IdealBranch, ...]
    particles: int
    degrees_of_freedom: int
    hot_temperature: float
    cold_temperature: float
    efficiency: float
    work_by_gas: float


def ideal_carnot_reference(simulation, points_per_branch: int = 65
                           ) -> IdealCarnotReference | None:
    """Return P=N*T/A curves; no measured pressure or temperature is fitted.

    The reference assumes point particles, thermodynamic equilibrium, and
    adiabats with T*A**(2/f) constant. Mixed rotational degrees of freedom or
    an area profile inconsistent with those adiabats has no such reference.
    """
    if points_per_branch < 2:
        raise ValueError("need at least two points per branch")
    mechanism = simulation.world.mechanism
    cam = getattr(mechanism, "cam", None)
    if cam is None:
        return None
    hot = getattr(mechanism, "hot_temperature", None)
    cold = getattr(mechanism, "cold_temperature", None)
    if hot is None or cold is None or not hot > cold > 0:
        return None
    inertia = simulation.world.bodies.inertia
    if np.all(np.isfinite(inertia)):
        dof = 3
    elif np.all(~np.isfinite(inertia)):
        dof = 2
    else:
        return None
    a1, a2, a3, a4 = cam.areas
    if min(cam.areas) <= 0:
        return None
    adiabatic_ratio = (hot/cold)**(dof/2)
    if not (math.isclose(a3/a2, adiabatic_ratio, rel_tol=1e-10) and
            math.isclose(a4/a1, adiabatic_ratio, rel_tol=1e-10)):
        return None
    n = simulation.world.bodies.n
    names = ("hot", "adiabatic_expansion", "cold", "adiabatic_compression")
    branches = []
    boundaries = cam.boundaries
    for index, name in enumerate(names):
        area = tuple(cam.area(float(boundaries[index]+(
            boundaries[index+1]-boundaries[index])*j/(points_per_branch-1)))
            for j in range(points_per_branch))
        if index == 0:
            temperature = (hot for _ in area)
        elif index == 1:
            temperature = (hot*(a2/a)**(2/dof) for a in area)
        elif index == 2:
            temperature = (cold for _ in area)
        else:
            temperature = (cold*(a4/a)**(2/dof) for a in area)
        pressure = tuple(n*t/a for t,a in zip(temperature,area))
        branches.append(IdealBranch(name,area,pressure))
    work = n*(hot-cold)*math.log(a2/a1)
    return IdealCarnotReference(tuple(branches),n,dof,hot,cold,
                                1-cold/hot,work)
