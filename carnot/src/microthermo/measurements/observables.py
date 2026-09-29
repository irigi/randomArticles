from __future__ import annotations

import numpy as np
from ..core.state import BodyArrays


def translational_temperature(state: BodyArrays, indices: np.ndarray | None = None) -> float:
    ids = np.arange(state.n) if indices is None else np.asarray(indices)
    if len(ids) < 2:
        return float("nan")
    masses = state.mass[ids]
    mean = np.sum(masses[:, None]*state.vel[ids], axis=0)/np.sum(masses)
    energy = 0.5*np.sum(masses[:, None]*(state.vel[ids]-mean)**2)
    dof = 2*len(ids)-2
    return float(2*energy/dof)


def rotational_temperature(state: BodyArrays, indices: np.ndarray | None = None) -> float:
    ids = np.arange(state.n) if indices is None else np.asarray(indices)
    ids = ids[np.isfinite(state.inertia[ids])]
    if not len(ids):
        return float("nan")
    return float(np.sum(state.inertia[ids]*state.omega[ids]**2)/len(ids))


def pressure(total_normal_impulse: float, length: float, elapsed: float) -> float:
    return total_normal_impulse/(length*elapsed) if length > 0 and elapsed > 0 else float("nan")


def chemical_potential_delta(concentration: float, reference: float) -> float:
    if concentration <= 0 or reference <= 0:
        return float("-inf")
    return float(np.log(concentration/reference))

