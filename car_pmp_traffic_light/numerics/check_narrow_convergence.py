"""Convergence check for the narrow-prior anticipatory-acceleration result."""
from model import Model
from solve_policies import solve_uniform
m=Model(900)
print('N, unrestricted, constrained, advantage')
for N in (275,550,825,1100):
    a=solve_uniform(27,27.5,N=N,model=m,start_labels=['legacy','brake_recover'])
    c=solve_uniform(27,27.5,N=N,model=m,no_speed_increasing_throttle=True,start_labels=['legacy','brake_recover'])
    print(N,a['cost'],c['cost'],c['cost']-a['cost'])
