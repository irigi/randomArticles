"""Scene bounds and aspect-preserving camera math, independent of Qt."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Bounds:
    left: float
    bottom: float
    right: float
    top: float


@dataclass(frozen=True)
class Camera:
    bounds: Bounds
    width: int
    height: int
    padding: int = 20

    @property
    def scale(self) -> float:
        span_x = max(self.bounds.right-self.bounds.left, 1e-9)
        span_y = max(self.bounds.top-self.bounds.bottom, 1e-9)
        return max(1e-9, min((self.width-2*self.padding)/span_x,
                             (self.height-2*self.padding)/span_y))

    def map(self, x: float, y: float) -> tuple[float, float]:
        center_x = (self.bounds.left+self.bounds.right)/2
        center_y = (self.bounds.bottom+self.bounds.top)/2
        return (self.width/2+(x-center_x)*self.scale,
                self.height/2-(y-center_y)*self.scale)


def scene_bounds(simulation) -> Bounds:
    """Use a fixed full-cycle envelope for Carnot and current geometry elsewhere."""
    mechanism = simulation.world.mechanism
    if hasattr(mechanism, "cam") and hasattr(mechanism, "apparatus_state"):
        cam = mechanism.cam
        width = max(cam.areas)/cam.height
        height = cam.height
        return Bounds(-.15*height, -.85*height,
                      width+.15*height, 1.15*height)

    coordinates = []
    for wall in simulation.world.walls:
        if hasattr(wall, "start"):
            coordinates.extend((tuple(wall.start),tuple(wall.end)))
        else:
            nx,ny = wall.inward_normal
            tx,ty = -ny,nx
            x,y = wall.point
            half = wall.length/2
            coordinates.extend(((float(x+half*tx),float(y+half*ty)),
                                (float(x-half*tx),float(y-half*ty))))
    state = simulation.world.bodies
    for (x,y),radius in zip(state.pos,state.radius):
        coordinates.extend(((float(x-radius),float(y-radius)),
                            (float(x+radius),float(y+radius))))
    if not coordinates:
        return Bounds(-1.,-1.,1.,1.)
    xs,ys = zip(*coordinates)
    left,right = min(xs),max(xs)
    bottom,top = min(ys),max(ys)
    pad = .06*max(right-left,top-bottom,1.)
    return Bounds(left-pad,bottom-pad,right+pad,top+pad)
