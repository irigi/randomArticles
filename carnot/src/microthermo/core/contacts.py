from __future__ import annotations

import math
from itertools import combinations
import numpy as np

from .events import Contact
from .geometry import cross2
from .state import BodyArrays


def contact_velocity(state: BodyArrays, i: int, r: np.ndarray) -> np.ndarray:
    return state.vel[i] + state.omega[i] * np.array([-r[1], r[0]])


def inverse_effective_mass(state: BodyArrays, contact: Contact) -> float:
    n, p = np.asarray(contact.normal), np.asarray(contact.point)
    total = 0.0
    for i in (contact.a, contact.b):
        if i is None or not state.dynamic[i]:
            continue
        r = p - state.pos[i]
        total += 1.0/state.mass[i]
        if math.isfinite(state.inertia[i]):
            total += cross2(r, n)**2/state.inertia[i]
    return total


def relative_normal_velocity(state: BodyArrays, contact: Contact) -> float:
    n, p = np.asarray(contact.normal), np.asarray(contact.point)
    va = contact_velocity(state, contact.a, p-state.pos[contact.a])
    vb = np.zeros(2) if contact.b is None else contact_velocity(state, contact.b, p-state.pos[contact.b])
    return float((vb-va) @ n)


def apply_impulse(state: BodyArrays, contact: Contact, impulse: float) -> None:
    n, p = np.asarray(contact.normal), np.asarray(contact.point)
    for i, sign in ((contact.a, -1.0), (contact.b, 1.0)):
        if i is None or not state.dynamic[i]:
            continue
        jv = sign * impulse * n
        state.vel[i] += jv/state.mass[i]
        if math.isfinite(state.inertia[i]):
            state.omega[i] += cross2(p-state.pos[i], jv)/state.inertia[i]
        state.version[i] += 1


def resolve_elastic(state: BodyArrays, contact: Contact, restitution: float = 1.0) -> tuple[float, float]:
    g = relative_normal_velocity(state, contact)
    if g >= 0.0:
        return 0.0, g
    d = inverse_effective_mass(state, contact)
    if d <= 0.0:
        raise ValueError("contact has no dynamic degree of freedom")
    impulse = -(1.0 + restitution)*g/d
    apply_impulse(state, contact, impulse)
    return impulse, g


def elastic_cluster_impulses(state: BodyArrays, contacts: tuple[Contact, ...],
                             velocity_tolerance: float = 1e-12) -> tuple[np.ndarray, float, int]:
    """Solve a small simultaneous stationary elastic cluster without changing state.

    Each contact's normal points from body a toward body b or a fixed wall.
    The returned nonnegative impulses reverse active normal velocities, leave
    all other contact velocities nonnegative, and conserve kinetic energy.
    """
    count = len(contacts)
    if count < 2 or count > 6:
        raise ValueError("elastic cluster size must be between 2 and 6")
    rows = []
    for contact in contacts:
        row = np.zeros((state.n, 3), dtype=float)
        n = np.asarray(contact.normal, float)
        p = np.asarray(contact.point, float)
        for i, sign in ((contact.a, -1.0), (contact.b, 1.0)):
            if i is None or not state.dynamic[i]:
                continue
            row[i, :2] += sign*n
            row[i, 2] += cross2(p-state.pos[i], sign*n)
        rows.append(row)
    jacobian = np.asarray(rows).reshape(count, 3*state.n)
    inverse_mass = np.empty((state.n, 3), dtype=float)
    inverse_mass[:, :2] = 1.0/state.mass[:, None]
    inverse_mass[:, 2] = np.where(np.isfinite(state.inertia), 1.0/state.inertia, 0.0)
    inverse_mass[~state.dynamic] = 0.0
    gram = (jacobian*inverse_mass.ravel()) @ jacobian.T
    velocity = np.concatenate((state.vel, state.omega[:, None]), axis=1).ravel()
    incoming = jacobian @ velocity
    scale = max(1.0, float(np.max(np.abs(incoming))))
    admissible = velocity_tolerance*scale
    for size in range(count, -1, -1):
        for active in combinations(range(count), size):
            impulses = np.zeros(count)
            rank = 0
            if active:
                matrix = gram[np.ix_(active, active)]
                rhs = -2.0*incoming[list(active)]
                solved, _, rank, _ = np.linalg.lstsq(matrix, rhs, rcond=1e-12)
                if np.max(np.abs(matrix@solved-rhs)) > admissible:
                    continue
                impulses[list(active)] = solved
            if np.min(impulses) < -admissible:
                continue
            impulses = np.maximum(impulses, 0.0)
            outgoing = incoming+gram@impulses
            if np.min(outgoing) < -admissible:
                continue
            energy_change = float(incoming@impulses + .5*impulses@gram@impulses)
            energy_scale = max(1.0, float(np.sum(np.abs(incoming*impulses))))
            if abs(energy_change) > 1e-10*energy_scale:
                continue
            residual = max(float(np.max(np.maximum(-outgoing, 0.0))),
                           abs(energy_change)/energy_scale)
            return impulses, residual, rank
    raise ValueError("no admissible energy-conserving elastic cluster solution")


def resolve_energy_step(state: BodyArrays, contact: Contact, delta_u: float,
                        threshold: float = 1e-12) -> tuple[bool, float]:
    g = relative_normal_velocity(state, contact)
    d = inverse_effective_mass(state, contact)
    radicand = g*g - 2.0*d*delta_u
    if radicand < -threshold:
        impulse = -2.0*g/d
        apply_impulse(state, contact, impulse)
        return False, impulse
    gp = math.copysign(math.sqrt(max(0.0, radicand)), g)
    impulse = (gp-g)/d
    apply_impulse(state, contact, impulse)
    return True, impulse
