from __future__ import annotations

import numpy as np
from .base import Experiment
from .common import box_walls, particle_specs
from ..config import RunConfig
from ..core.boundaries import BoundaryKind, SegmentWall
from ..core.state import BodyArrays, BodySpec, Shape
from ..core.geometry import closest_point_segment
from ..runner.simulation import World


class PermeabilityExperiment(Experiment):
    def build(self, config: RunConfig, rng: np.random.Generator) -> World:
        width,height=4.,2.
        small=particle_specs(config.particles,width,height,.035,config.temperature,rng)
        hosts=[]
        for y in (.35,.75,1.15,1.65):
            hosts.append(BodySpec((3.25,y),tuple(rng.normal(0,.18,2)),mass=20.,radius=.16,
                                  shape=Shape.DISC,species=1))
        walls=box_walls(width,height,BoundaryKind.HOT,config.temperature)
        x=2.; gap=.11; pitch=.28
        cursor=0.
        idx=0
        while cursor<height:
            end=min(height,cursor+pitch-gap)
            if end>cursor: walls.append(SegmentWall((x,cursor),(x,end),.018,name=f"pore_wall_{idx}"))
            cursor+=pitch; idx+=1
        segments=[w for w in walls if isinstance(w,SegmentWall)]
        for i,spec in enumerate(small):
            for _ in range(10000):
                p=np.asarray(spec.position)
                clear=all(np.linalg.norm(p-closest_point_segment(p,w.start,w.end)[0])>.035+w.thickness+1e-8
                          for w in segments)
                clear &= all(np.linalg.norm(p-np.asarray(host.position)) > spec.radius+host.radius+1e-8
                             for host in hosts)
                clear &= all(j == i or np.linalg.norm(p-np.asarray(other.position)) > spec.radius+other.radius+1e-8
                             for j,other in enumerate(small))
                if clear: break
                spec.position=(float(rng.uniform(.035,width-.035)),float(rng.uniform(.035,height-.035)))
            else:
                raise ValueError("could not place discs without overlap")
        return World(BodyArrays.from_specs(small+hosts),walls,
                     metadata={"name":"selective_membrane","gap":gap,
                               "disc_diameter":.07,"host_diameter":.32,
                               "interpretation":"geometric selective permeability"})
