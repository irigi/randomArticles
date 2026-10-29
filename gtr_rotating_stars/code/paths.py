"""Shared project paths, so scripts work from any working directory."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'            # generated results
INPUTS = DATA / 'inputs'        # fixed reference inputs (not regenerated)
FIG = ROOT / 'figures'          # figures included by the LaTeX sources
DATA.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)
