"""Regenerate the retained rotation and stability figures from baseline scripts."""
from pathlib import Path
import subprocess,sys,shutil
root=Path(__file__).resolve().parent
subprocess.run([sys.executable,'rewrite_figures.py'],cwd=root,check=True)
for name in ('rotation','stability'):
    for ext in ('.pdf','.png'):
        shutil.copy2(root/'output/pdf/revision/figures'/(name+ext),root/'output/pdf/final/figures'/(name+ext))
