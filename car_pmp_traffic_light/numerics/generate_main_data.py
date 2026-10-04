"""Regenerate numerical policy arrays used by the paper.
Run from this directory or any directory; outputs go to ../data.
"""
from pathlib import Path
import numpy as np
from model import Model, published_accel_cap
from solve_policies import solve_known, solve_uniform, save_npz

HERE=Path(__file__).resolve().parent
OUT=HERE.parent/'data'; OUT.mkdir(exist_ok=True)

base=Model(300)
ALL_UNCERTAIN=['legacy','rolling_0.65','rolling_1.05','stop_wait','coast_then_stop','brake_recover']
BROAD_BOUNDED=['legacy','rolling_1.05','stop_wait','coast_then_stop','brake_recover']
NARROW_BOUNDED=['legacy','brake_recover']
for name,r in [
    ('k_28.npz', solve_known(28,N=336,model=base)),
    ('u_0_60.npz', solve_uniform(0,60,N=480,model=base,start_labels=BROAD_BOUNDED)),
    ('u_26_30.npz', solve_uniform(26,30,N=360,model=base,start_labels=ALL_UNCERTAIN)),
    ('narrow300_N1100.npz', solve_uniform(27,27.5,N=1100,model=base,start_labels=NARROW_BOUNDED)),
    ('narrow300_noinc_N1100.npz', solve_uniform(27,27.5,N=1100,model=base,no_speed_increasing_throttle=True,start_labels=NARROW_BOUNDED)),
]: save_npz(OUT/name,r)

m900=Model(900)
for name,r in [
    ('narrow900_N1100.npz', solve_uniform(27,27.5,N=1100,model=m900,start_labels=NARROW_BOUNDED)),
    ('narrow900_noacc_N1100.npz', solve_uniform(27,27.5,N=1100,model=m900,no_speed_increasing_throttle=True,start_labels=NARROW_BOUNDED)),
]: save_npz(OUT/name,r)

# Mechanical sensitivity: acceleration cap calibrated to 8.4 s 0-100 km/h,
# and brake addition reduced from 7 to 4 m/s^2.
cap=published_accel_cap()
ms=Model(300,accel_cap=cap,brake_a=4.0)
for name,r in [
    ('sens_known28.npz',solve_known(28,N=336,model=ms)),
    ('sens_broad60.npz',solve_uniform(0,60,N=480,model=ms,start_labels=BROAD_BOUNDED)),
    ('sens_late.npz',solve_uniform(26,30,N=360,model=ms,start_labels=ALL_UNCERTAIN)),
]: save_npz(OUT/name,r)
ms900=Model(900,accel_cap=cap,brake_a=4.0)
for name,r in [
    ('sens_narrow900.npz',solve_uniform(27,27.5,N=550,model=ms900,start_labels=NARROW_BOUNDED)),
    ('sens_narrow900_noinc.npz',solve_uniform(27,27.5,N=550,model=ms900,no_speed_increasing_throttle=True,start_labels=NARROW_BOUNDED)),
]: save_npz(OUT/name,r)
print('wrote',OUT)
