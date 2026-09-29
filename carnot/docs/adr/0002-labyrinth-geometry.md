# ADR 0002: finite segment labyrinth reference geometry

Status: accepted with extension required, 2026-09-28.

Labyrinth layers use finite two-sided capsule-like segments with explicit wall
thickness and positive aperture clearance. Oriented portal surfaces are
separate event geometry. Their crossings update region membership and, when
configured, apply reversible potential-step impulses.

This representation provides real openings, passage-edge impacts, and hard
capture without treating union seams as physical. The initial reference preset
fixes these segments to the laboratory support. A production mobile host needs
one shared rigid state for all segments, exposed-union filtering, and a per-host
BVH before it can be used for the specification's recoil and torque studies.

