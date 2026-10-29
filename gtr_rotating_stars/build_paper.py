"""Assemble editable components and compile PDFs with a normal TeX installation."""
from pathlib import Path
import subprocess,argparse
parser=argparse.ArgumentParser()
parser.add_argument('--assemble-only',action='store_true')
args=parser.parse_args()
root=Path(__file__).resolve().parent
paper=root
body=(paper/'article_body.tex').read_text()
bib=(paper/'references.tex').read_text()
abstract,main=body.split(r'\end{abstract}',1)
for stem in ('journal','companion'):
    prefix=(paper/f'{stem}_preamble.tex').read_text()
    if stem=='journal':
        text=prefix+abstract+r'\end{abstract}'+'\n'+r'\maketitle'+'\n'+main+'\n'+r'\clearpage'+'\n'+bib+'\n'+r'\end{document}'+'\n'
    else:
        text=prefix+body+'\n'+bib+'\n'+r'\end{document}'+'\n'
    (paper/f'rotating_star_{stem}.tex').write_text(text)
if not args.assemble_only:
    for stem in ('journal','companion','supplement'):
        for _ in range(3):
            subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error',f'rotating_star_{stem}.tex'],cwd=paper,check=True)
