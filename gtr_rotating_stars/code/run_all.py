"""Regenerate all data in data/ and all figures in figures/, in dependency order."""
import subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = [
    'technical_family.py',              # generating family, thresholds
    'technical_stability.py',           # fixed-EOS sequences
    'technical_radial.py',              # radial modes
    'technical_verification.py',        # Einstein-equation residuals
    'rewrite_comparison.py',            # reference-model comparison
    'benchmark.py',                     # external benchmarks
    'thermodynamics.py',                # thermodynamic profiles
    'figures.py',                       # seven article/supplement figures
    'figures_rotation_stability.py',    # rotation and stability figures
]
for step in STEPS:
    t = time.time()
    print(f'== {step}', flush=True)
    subprocess.run([sys.executable, str(HERE / step)], cwd=HERE, check=True)
    print(f'   done in {time.time() - t:.0f} s', flush=True)
