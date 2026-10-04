from __future__ import annotations

from dataclasses import dataclass, asdict
from enum import Enum
from typing import Any


class TOIStatus(str, Enum):
    COLLISION = "collision"
    BRANCH = "branch"
    NO_COLLISION = "no_collision"
    INDETERMINATE = "indeterminate"


@dataclass(frozen=True)
class Contact:
    a: int
    b: int | None
    point: tuple[float, float]
    normal: tuple[float, float]
    feature_a: int = -1
    feature_b: int = -1
    boundary: int | None = None


@dataclass(frozen=True)
class TOIResult:
    status: TOIStatus
    time: float | None = None
    error: float = 0.0
    contact: Contact | None = None
    reason: str = ""
    contacts: tuple[Contact, ...] = ()
    transition: dict[str, Any] | None = None


@dataclass
class InteractionRecord:
    time: float
    kind: str
    participants: tuple[int, ...]
    impulse: tuple[float, float] = (0.0, 0.0)
    energy_before: float = 0.0
    energy_after: float = 0.0
    heat_into_system: float = 0.0
    work_on_system: float = 0.0
    metadata: dict[str, Any] | None = None
    contact_point: tuple[float, float] | None = None

    def as_dict(self) -> dict[str, Any]:
        record = asdict(self)
        if self.contact_point is None:
            record.pop("contact_point")
        return record
