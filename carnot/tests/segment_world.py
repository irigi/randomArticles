"""Segment-wall partition world for reference-engine tests.

This is the geometry of the former ``selective_membrane`` stub preset: four
heavy discs left of a partition of finite segment walls, with small discs
anywhere clear of them. It exercises the reference engine's segment-wall
paths, which the osmosis presets (edmd posts) no longer use.
"""
import numpy as np

from microthermo.core.boundaries import BoundaryKind, SegmentWall
from microthermo.core.geometry import closest_point_segment
from microthermo.core.state import BodyArrays, BodySpec, Shape
from microthermo.experiments.common import box_walls, particle_specs
from microthermo.runner.simulation import World


def segment_partition_world(particles: int, seed: int = 123) -> World:
    rng = np.random.default_rng(seed)
    width, height = 4., 2.
    small = particle_specs(particles, width, height, .035, 1., rng)
    hosts = [BodySpec((3.25, y), tuple(rng.normal(0, .18, 2)), mass=20., radius=.16,
                      shape=Shape.DISC, species=1) for y in (.35, .75, 1.15, 1.65)]
    walls = box_walls(width, height, BoundaryKind.HOT, 1.)
    x, gap, pitch, cursor, index = 2., .11, .28, 0., 0
    while cursor < height:
        end = min(height, cursor + pitch - gap)
        if end > cursor:
            walls.append(SegmentWall((x, cursor), (x, end), .018, name=f"pore_wall_{index}"))
        cursor += pitch
        index += 1
    segments = [w for w in walls if isinstance(w, SegmentWall)]
    for i, spec in enumerate(small):
        for _ in range(10000):
            p = np.asarray(spec.position)
            clear = all(np.linalg.norm(p - closest_point_segment(p, w.start, w.end)[0])
                        > .035 + w.thickness + 1e-8 for w in segments)
            clear &= all(np.linalg.norm(p - np.asarray(h.position)) > spec.radius + h.radius + 1e-8
                         for h in hosts)
            clear &= all(j == i or np.linalg.norm(p - np.asarray(o.position))
                         > spec.radius + o.radius + 1e-8 for j, o in enumerate(small))
            if clear:
                break
            spec.position = (float(rng.uniform(.035, width - .035)),
                             float(rng.uniform(.035, height - .035)))
        else:
            raise ValueError("could not place discs without overlap")
    return World(BodyArrays.from_specs(small + hosts), walls,
                 metadata={"name": "segment_partition"})
