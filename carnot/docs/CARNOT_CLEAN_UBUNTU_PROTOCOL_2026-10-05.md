# Clean Ubuntu release check — prepared 2026-10-05

Section 5 of `CARNOT_GAP_CLOSURE_PLAN.md` still needs installation evidence
from a clean Ubuntu image. The earlier fresh virtual environment passed on an
existing Ubuntu 24.04.5 host; it does not establish what system packages a
clean image may need. Run `bash runs/release_clean_ubuntu.sh` from a
checkout of this repository on a clean Ubuntu 24.04 system with Python 3.12,
`python3-venv`, Git, and network access to package wheels. The launcher does
not install system packages or require a display.

The launcher records OS, Python, kernel, CPU, Git revision, and repository
status. In parallel, it tests two isolated installations from committed source
archives:

1. Headless package with `requirements-lock.txt`, `pip install .`, a check
   that importing `microthermo` does not load Qt, and quick validation.
2. GUI plus accelerator with `requirements-accel-lock.txt`, `pip install
   '.[gui,accel]'`, the complete unittest suite on the offscreen Qt platform,
   and scientific validation.

Each job saves stdout, stderr, exit code, and installed package versions under
its own name. The batch directory is temporary; after review, retain a compact
installation report and delete the batch. A pass on a clean image
supports installation and offscreen GUI compatibility. The previous real
display acceptance applies to the measured xcb host, not automatically to
another display backend or GPU.

The first run of this launcher on 2026-10-05 used the **existing** ThinkPad
host, not a clean OS image. Both isolated installs passed: the headless import
loaded no Qt and all 16 quick checks passed; the GUI/accelerator install
passed 157 unit tests and all 16 scientific checks. The compact environment
and package summary is
`runs/gap_followup_2026-10-04/release_archive_install_2026-10-05.json`.
This validates independent installs from committed source but leaves the
clean-image gate open. The launcher is retained as `runs/release_clean_ubuntu.sh`
for use on an actual clean Ubuntu machine.
