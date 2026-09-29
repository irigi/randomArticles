from __future__ import annotations

import numpy as np
from .base import Experiment
from .common import box_walls, particle_specs
from ..config import RunConfig
from ..core.boundaries import BoundaryKind, SegmentWall, Wall
from ..core.state import BodyArrays
from ..core.geometry import closest_point_segment
from ..runner.simulation import Portal, World


class LabyrinthExperiment(Experiment):
    """Hard-wall nested chambers whose occupancy follows trajectory crossings."""

    def __init__(self, binding_energy: float = 0.0):
        self.binding_energy=binding_energy

    def build(self, config: RunConfig, rng: np.random.Generator) -> World:
        width, height, radius = 4.0, 2.0, 0.035
        specs = particle_specs(config.particles, width, height, radius,
                               config.temperature, rng)
        walls = box_walls(width, height, BoundaryKind.HOT, config.temperature)
        portals: list[Portal] = []
        # Three nested, offset C-shaped rectangular chambers.  Segments are
        # finite in the renderer/export; the reference collision walls use
        # their supporting lines and therefore this preset uses one-sided
        # chamber surfaces with generous clearances.
        centers = [(2.7, 1.0), (2.7, 1.0), (2.68, 1.03)]
        sizes = [(0.8, 0.7), (0.55, 0.45), (0.3, 0.23)]
        gaps = [("left", 0.25), ("top", 0.18), ("right", 0.12)]
        for layer, ((cx, cy), (w, h), (side, gap)) in enumerate(zip(centers, sizes, gaps), 1):
            # Portal surfaces are exact occupancy measurements. Geometry is
            # included as metadata for clients with finite-segment collision kernels.
            if side == "left": p, n = np.array([cx-w/2, cy]), np.array([-1., 0.])
            elif side == "right": p, n = np.array([cx+w/2, cy]), np.array([1., 0.])
            else: p, n = np.array([cx, cy+h/2]), np.array([0., 1.])
            level=-self.binding_energy if layer==3 else 0.0
            outer_level=0.0
            portals.append(Portal(p, n, gap/2, host_id=0, inside_label=layer,
                                  outside_label=layer-1, delta_u=level-outer_level,
                                  inside_energy=level))
            left,right,bottom,top=cx-w/2,cx+w/2,cy-h/2,cy+h/2
            def split_segment(start,end,opening_axis,opening_center,opening_width):
                lo,hi=opening_center-opening_width/2,opening_center+opening_width/2
                if opening_axis==0:
                    return [(start,(lo,start[1])),((hi,start[1]),end)]
                return [(start,(start[0],lo)),((start[0],hi),end)]
            segments=[((left,bottom),(right,bottom)),((left,top),(right,top)),
                      ((left,bottom),(left,top)),((right,bottom),(right,top))]
            if side=="left": segments=segments[:2]+split_segment((left,bottom),(left,top),1,cy,gap)+segments[3:]
            elif side=="right": segments=segments[:3]+split_segment((right,bottom),(right,top),1,cy,gap)
            elif side=="top": segments=[segments[0]]+split_segment((left,top),(right,top),0,cx,gap)+segments[2:]
            for j,(start,end) in enumerate(segments):
                if np.linalg.norm(np.asarray(end)-start)>1e-12:
                    walls.append(SegmentWall(start,end,.012,name=f"host0_layer{layer}_wall{j}"))
        # Initial conditions may occupy chambers, but never material walls.
        segment_walls=[w for w in walls if isinstance(w,SegmentWall)]
        for i,spec in enumerate(specs):
            for _ in range(10000):
                p=np.asarray(spec.position)
                clear=all(np.linalg.norm(p-closest_point_segment(p,w.start,w.end)[0]) > radius+w.thickness+1e-8
                          for w in segment_walls)
                separate=all(j==i or np.linalg.norm(p-np.asarray(other.position))>2.05*radius
                             for j,other in enumerate(specs))
                if clear and separate: break
                spec.position=(float(rng.uniform(radius,width-radius)),float(rng.uniform(radius,height-radius)))
            else: raise ValueError("could not create overlap-free labyrinth state")
        name="labyrinth_energetic" if self.binding_energy else "labyrinth"
        return World(BodyArrays.from_specs(specs), walls, portals=portals,
                     metadata={"name": name, "binding_energy":self.binding_energy,
                        "chambers": [dict(center=c, size=s, gap=g)
                        for c, s, g in zip(centers, sizes, gaps)],
                        "interpretation": "hard-only geometric association",
                        "accessible_pocket_areas": [w*h for w, h in sizes]})

    def instruments(self):
        return super().instruments()+["occupancy", "entry_exit", "residence_time", "chemical_potential_reference"]
