from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
import math
import numpy as np

_GAUSS_NODES, _GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(8)


def smooth5(u: float) -> float:
    u = min(1.0, max(0.0, u))
    return u**3 * (10.0 + u*(-15.0 + 6.0*u))


def smooth5_derivative(u: float) -> float:
    u = min(1.0, max(0.0, u))
    return 30.0*u*u*(1.0-u)*(1.0-u)


@dataclass(frozen=True)
class CarnotCam:
    """C2 periodic area profile and right-continuous thermal selector."""

    areas: tuple[float, float, float, float]
    fractions: tuple[float, float, float, float] = (0.25, 0.25, 0.25, 0.25)
    height: float = 1.0

    @classmethod
    def design(cls, a1: float, ratio: float, t_hot: float, t_cold: float,
               degrees_of_freedom: int, height: float = 1.0) -> "CarnotCam":
        if not (a1 > 0 and ratio > 1 and t_hot > t_cold > 0):
            raise ValueError("require A1>0, ratio>1, and T_hot>T_cold>0")
        alpha = (t_hot/t_cold)**(degrees_of_freedom/2)
        return cls((a1, a1*ratio, a1*ratio*alpha, a1*alpha), height=height)

    @cached_property
    def boundaries(self) -> np.ndarray:
        boundaries = np.concatenate(([0.0], np.cumsum(self.fractions))) * (2*np.pi)
        boundaries.flags.writeable = False
        return boundaries

    def _sector(self, phi: float) -> tuple[int, float, float]:
        phase = phi % (2*np.pi)
        b = self.boundaries
        i = min(3, int(np.searchsorted(b, phase, side="right")-1))
        width = b[i+1]-b[i]
        return i, (phase-b[i])/width, width

    def area(self, phi: float) -> float:
        i, u, _ = self._sector(phi)
        return self.areas[i] + (self.areas[(i+1)%4]-self.areas[i])*smooth5(u)

    def area_derivative(self, phi: float) -> float:
        i, u, width = self._sector(phi)
        return (self.areas[(i+1)%4]-self.areas[i])*smooth5_derivative(u)/width

    def piston_x(self, phi: float) -> float:
        return self.area(phi)/self.height

    def piston_dx_dphi(self, phi: float) -> float:
        return self.area_derivative(phi)/self.height

    def branch(self, phi: float, reversed_cycle: bool = False) -> str:
        names = ("hot", "adiabatic_expansion", "cold", "adiabatic_compression")
        i = self._sector(-phi if reversed_cycle else phi)[0]
        return names[i]


@dataclass
class Shaft:
    phi: float
    momentum: float
    inertia: float
    prescribed_omega: float | None = None
    load_torque: float = 0.0
    spring_k: float = 0.0
    spring_rest: float = 0.0
    piston_mass: float = 0.0
    cam: CarnotCam | None = None

    def __post_init__(self) -> None:
        if self.inertia <= 0 or self.spring_k < 0 or self.load_torque < 0 or self.piston_mass < 0:
            raise ValueError("require positive shaft inertia and nonnegative masses, spring, and load")
        if self.piston_mass and self.cam is None:
            raise ValueError("reflected piston mass requires a cam")

    def effective_inertia(self, phi: float | None = None) -> float:
        phase = self.phi if phi is None else phi
        if not self.piston_mass:
            return self.inertia
        slope = self.cam.piston_dx_dphi(phase)
        return self.inertia+self.piston_mass*slope*slope

    def _metric_distance(self, start: float, end: float) -> float:
        """Signed arc length for the one-coordinate kinetic metric."""
        width = end-start
        pieces = max(1,math.ceil(abs(width)/.1))
        total = 0.0
        for j in range(pieces):
            left = start+width*j/pieces
            right = start+width*(j+1)/pieces
            middle,half = .5*(left+right),.5*(right-left)
            total += half*sum(float(weight)*math.sqrt(self.effective_inertia(
                middle+half*float(node))) for node,weight in
                zip(_GAUSS_NODES,_GAUSS_WEIGHTS))
        return total

    def _free_reflected_trajectory(self, dt: float) -> tuple[float, float, float]:
        mass0 = self.effective_inertia()
        metric_speed = self.momentum/math.sqrt(mass0)
        if metric_speed == 0 or dt == 0:
            return self.phi,self.momentum,0.0
        direction = math.copysign(1.0,metric_speed)
        reach = abs(metric_speed)*dt/math.sqrt(self.inertia)
        target = abs(metric_speed)*dt
        lo,hi = 0.0,reach
        for _ in range(60):
            middle = .5*(lo+hi)
            distance = abs(self._metric_distance(self.phi,self.phi+direction*middle))
            if distance >= target:
                hi = middle
            else:
                lo = middle
            if hi-lo <= 2e-14*max(1.0,reach):
                break
        phi = self.phi+direction*.5*(lo+hi)
        mass = self.effective_inertia(phi)
        return phi,metric_speed*math.sqrt(mass),metric_speed/math.sqrt(mass)

    def _coupled_step(self, phi: float, momentum: float, dt: float
                      ) -> tuple[float, float, float, float | None]:
        """Scalar discrete-gradient step, split at a load/spring turning point."""
        mass0 = self.effective_inertia(phi)
        spring_force = -self.spring_k*(phi-self.spring_rest)
        if abs(momentum) <= 1e-14:
            if abs(spring_force) <= self.load_torque:
                return phi,0.0,0.0,None
            direction = math.copysign(1.0,spring_force)
            momentum = 0.0
        else:
            direction = math.copysign(1.0,momentum)
        energy = (momentum*momentum/(2*mass0)+
                  .5*self.spring_k*(phi-self.spring_rest)**2)

        def available(travel: float) -> float:
            next_phi = phi+direction*travel
            return energy-self.load_torque*travel-.5*self.spring_k*(
                next_phi-self.spring_rest)**2

        upper = max(1e-8,2*dt*math.sqrt(2*energy/self.inertia))
        for _ in range(100):
            if available(upper) <= 0:
                break
            upper *= 2
        else:
            raise RuntimeError("could not bracket shaft turning angle")
        lower = 0.0
        for _ in range(65):
            middle = .5*(lower+upper)
            if available(middle) > 0:
                lower = middle
            else:
                upper = middle
        turn_travel = .5*(lower+upper)
        turn_mass = self.effective_inertia(phi+direction*turn_travel)
        if momentum:
            turn_speed = .25*abs(momentum)*(1/mass0+1/turn_mass)
            turn_time = turn_travel/turn_speed
        else:
            turn_time = math.inf
        if dt >= turn_time:
            return phi+direction*turn_travel,0.0,turn_travel,turn_time

        lower,upper = 0.0,turn_travel
        for _ in range(65):
            travel = .5*(lower+upper)
            next_phi = phi+direction*travel
            next_mass = self.effective_inertia(next_phi)
            next_momentum = math.sqrt(max(0.0,2*next_mass*available(travel)))
            mean_speed = .25*(abs(momentum)+next_momentum)*(
                1/mass0+1/next_mass)
            if travel > dt*mean_speed:
                upper = travel
            else:
                lower = travel
        travel = .5*(lower+upper)
        next_phi = phi+direction*travel
        next_mass = self.effective_inertia(next_phi)
        next_momentum = direction*math.sqrt(max(0.0,2*next_mass*available(travel)))
        return next_phi,next_momentum,travel,None

    def _coupled_trajectory(self, dt: float
                            ) -> tuple[float, float, float, tuple[float, ...]]:
        phi,momentum = self.phi,self.momentum
        elapsed = travel = 0.0
        turns: list[float] = []
        max_step = min(.01,.25*math.sqrt(self.inertia/self.spring_k)) if self.spring_k else .01
        for _ in range(10000):
            remaining = dt-elapsed
            if remaining <= 1e-15*max(1.0,dt):
                return phi,momentum,travel,tuple(turns)
            step = min(remaining,max_step)
            phi,momentum,segment,turn = self._coupled_step(phi,momentum,step)
            travel += segment
            elapsed += turn if turn is not None else step
            if turn is not None and elapsed < dt:
                turns.append(elapsed)
        raise RuntimeError("coupled shaft trajectory exceeded step limit")

    @property
    def omega(self) -> float:
        return (self.prescribed_omega if self.prescribed_omega is not None else
                self.momentum/self.effective_inertia())

    def kinetic_energy(self) -> float:
        return 0.5*self.effective_inertia()*self.omega**2

    def spring_energy(self) -> float:
        return 0.5*self.spring_k*(self.phi-self.spring_rest)**2

    def _loaded_trajectory(self, dt: float) -> tuple[float, float, float, tuple[float, ...]]:
        """Exact spring motion between Coulomb-load turns, with load travel."""
        phi, speed = self.phi, self.momentum/self.inertia
        elapsed = travel = 0.0
        turns: list[float] = []
        for _ in range(10000):
            remaining = dt-elapsed
            if remaining <= 0:
                return phi,self.inertia*speed,travel,tuple(turns)
            spring_force = -self.spring_k*(phi-self.spring_rest)
            if abs(speed) <= 1e-14:
                speed = 0.0
                if abs(spring_force) <= self.load_torque:
                    return phi,0.0,travel,tuple(turns)
                direction = math.copysign(1.0,spring_force)
            else:
                direction = math.copysign(1.0,speed)
            if self.spring_k == 0:
                acceleration = -direction*self.load_torque/self.inertia
                turn = -speed/acceleration if acceleration else math.inf
                step = min(remaining,turn)
                next_phi = phi+speed*step+0.5*acceleration*step*step
                next_speed = speed+acceleration*step
            else:
                frequency = math.sqrt(self.spring_k/self.inertia)
                center = self.spring_rest-direction*self.load_torque/self.spring_k
                displacement = phi-center
                phase = math.atan2(displacement*frequency,speed)
                turn_angle = (math.pi/2-phase) % math.pi
                if turn_angle <= 1e-14:
                    turn_angle = math.pi
                turn = turn_angle/frequency
                step = min(remaining,turn)
                cosine,sine = math.cos(frequency*step),math.sin(frequency*step)
                next_phi = center+displacement*cosine+speed*sine/frequency
                next_speed = -displacement*frequency*sine+speed*cosine
            travel += abs(next_phi-phi)
            phi,speed = next_phi,next_speed
            elapsed += step
            if step < turn or elapsed >= dt:
                return phi,self.inertia*speed,travel,tuple(turns)
            turns.append(elapsed)
            speed = 0.0
        raise RuntimeError("shaft load trajectory exceeded turning-point limit")

    def trajectory(self, dt: float) -> tuple[float, float, float]:
        """Position, momentum, and drift speed for the next free step.

        Collision queries and committed integration use this same map. The
        drift speed is the derivative of the predicted position with respect
        to the step length, rather than a constant-speed extrapolation.
        """
        if dt < 0:
            raise ValueError("shaft prediction requires a nonnegative interval")
        if self.prescribed_omega is not None:
            return self.phi+self.prescribed_omega*dt, self.momentum, self.prescribed_omega
        if self.piston_mass and (self.spring_k or self.load_torque):
            phi,momentum,_,_ = self._coupled_trajectory(dt)
            return phi,momentum,momentum/self.effective_inertia(phi)
        if self.piston_mass:
            return self._free_reflected_trajectory(dt)
        omega0 = self.momentum/self.inertia
        if self.load_torque == 0 and self.spring_k > 0:
            frequency = math.sqrt(self.spring_k/self.inertia)
            displacement = self.phi-self.spring_rest
            cosine, sine = math.cos(frequency*dt), math.sin(frequency*dt)
            phi = self.spring_rest+displacement*cosine+omega0*sine/frequency
            speed = -displacement*frequency*sine+omega0*cosine
            return phi, self.inertia*speed, speed
        if self.load_torque > 0:
            phi,momentum,_,_ = self._loaded_trajectory(dt)
            return phi,momentum,momentum/self.inertia
        torque0 = -self.spring_k*(self.phi-self.spring_rest)
        if omega0:
            torque0 -= math.copysign(self.load_torque, omega0)
        half_momentum = self.momentum+0.5*dt*torque0
        phi = self.phi+dt*half_momentum/self.inertia
        drift_speed = (self.momentum+dt*torque0)/self.inertia
        torque1 = -self.spring_k*(phi-self.spring_rest)
        if half_momentum:
            torque1 -= math.copysign(self.load_torque, half_momentum)
        momentum = half_momentum+0.5*dt*torque1
        return phi, momentum, drift_speed

    def speed_bound(self, dt: float) -> float:
        """Upper bound for the speed of the position predictor on [0, dt]."""
        if self.prescribed_omega is not None:
            return abs(self.prescribed_omega)
        if self.piston_mass and (self.spring_k or self.load_torque):
            energy = self.kinetic_energy()+self.spring_energy()
            return math.sqrt(2*energy/self.inertia)
        if self.piston_mass:
            metric_speed = self.momentum/math.sqrt(self.effective_inertia())
            return abs(metric_speed)/math.sqrt(self.inertia)
        omega0 = self.momentum/self.inertia
        if self.load_torque == 0 and self.spring_k > 0:
            return math.hypot(omega0, (self.phi-self.spring_rest)*
                              math.sqrt(self.spring_k/self.inertia))
        if self.load_torque > 0:
            return math.hypot(omega0, (self.phi-self.spring_rest)*
                              math.sqrt(self.spring_k/self.inertia))
        acceleration = (-self.spring_k*(self.phi-self.spring_rest)
                        - math.copysign(self.load_torque, omega0) if omega0 else
                        -self.spring_k*(self.phi-self.spring_rest))/self.inertia
        return max(abs(omega0), abs(omega0+acceleration*dt))

    def turning_times(self, dt: float) -> tuple[float, ...]:
        """Interior stationary times of the position predictor."""
        if self.prescribed_omega is not None or dt <= 0:
            return ()
        if self.piston_mass and (self.spring_k or self.load_torque):
            return self._coupled_trajectory(dt)[3]
        omega0 = self.momentum/self.inertia
        if self.load_torque == 0 and self.spring_k > 0:
            frequency = math.sqrt(self.spring_k/self.inertia)
            phase = math.atan2((self.phi-self.spring_rest)*frequency, omega0)
            first = (math.pi/2-phase)/frequency
            period = math.pi/frequency
            k0 = math.floor(-first/period)+1
            return tuple(t for k in range(k0, k0+math.ceil(dt/period)+2)
                         if 0 < (t := first+k*period) < dt)
        if self.load_torque > 0:
            return self._loaded_trajectory(dt)[3]
        acceleration = (-self.spring_k*(self.phi-self.spring_rest)
                        - math.copysign(self.load_torque, omega0) if omega0 else
                        -self.spring_k*(self.phi-self.spring_rest))/self.inertia
        if acceleration:
            turn = -omega0/acceleration
            if 0 < turn < dt:
                return (turn,)
        return ()

    def advance(self, dt: float) -> tuple[float, float]:
        old_phi = self.phi
        self.phi, self.momentum, _ = self.trajectory(dt)
        return old_phi, self.phi

    def advance_with_travel(self, dt: float) -> tuple[float, float, float]:
        old_phi = self.phi
        if (self.prescribed_omega is None and self.piston_mass and
            (self.spring_k or self.load_torque)):
            self.phi,self.momentum,travel,_ = self._coupled_trajectory(dt)
        elif self.prescribed_omega is None and self.load_torque > 0:
            self.phi,self.momentum,travel,_ = self._loaded_trajectory(dt)
        else:
            self.phi,self.momentum,_ = self.trajectory(dt)
            travel = abs(self.phi-old_phi)
        return old_phi,self.phi,travel


@dataclass
class LinearSpring:
    coordinate: float
    velocity: float
    mass: float
    stiffness: float
    rest: float = 0.0

    def energy(self) -> float:
        return 0.5*self.mass*self.velocity**2 + 0.5*self.stiffness*(self.coordinate-self.rest)**2

    def advance_exact(self, dt: float) -> None:
        w = math.sqrt(self.stiffness/self.mass)
        q = self.coordinate-self.rest
        c, s = math.cos(w*dt), math.sin(w*dt)
        self.coordinate = self.rest + q*c + self.velocity*s/w
        self.velocity = -q*w*s + self.velocity*c
