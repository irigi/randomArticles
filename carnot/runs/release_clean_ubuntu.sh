#!/usr/bin/env bash
# Run on a clean Ubuntu 24.04 checkout. Delete each generated
# TEMP_gap_release_clean_ubuntu_batch_* directory after review.
set -u
set -o pipefail
cd "$(dirname "$0")/.." || exit 1

batch=$(mktemp -d runs/TEMP_gap_release_clean_ubuntu_batch_XXXXXXXX)
printf '%s\n' 'TEMPORARY: delete after review and condensation.' > "$batch/TEMPORARY_DELETE_AFTER_REVIEW.txt"
{
  cat /etc/os-release
  python3 --version
  uname -a
  git rev-parse HEAD
  git status --short
  lscpu
} > "$batch/environment.stdout" 2> "$batch/environment.stderr"
printf 'Batch output: %s\n' "$batch"

if ! python3 -c 'import platform, sys, venv; o=platform.freedesktop_os_release(); assert o.get("ID") == "ubuntu" and o.get("VERSION_ID") == "24.04" and sys.version_info[:2] == (3, 12)' \
    > "$batch/prerequisites.stdout" 2> "$batch/prerequisites.stderr"; then
  printf '%s\n' 2 > "$batch/prerequisites.exit"
  printf 'Ubuntu 24.04 and Python 3.12 with venv are required; see %s/prerequisites.stderr\n' "$batch" >&2
  exit 2
fi
printf '%s\n' 0 > "$batch/prerequisites.exit"

run_job() {
  local label=$1 code
  mkdir -p "$batch/$label/source"
  if {
    git archive HEAD | tar -x -C "$batch/$label/source" &&
    python3 -m venv "$batch/$label/venv" &&
    (
      cd "$batch/$label/source" || exit 1
      if [[ "$label" == headless ]]; then
        ../venv/bin/python -m pip install --disable-pip-version-check \
          -r requirements-lock.txt &&
        ../venv/bin/python -m pip install --disable-pip-version-check . &&
        ../venv/bin/python -c \
          'import sys, microthermo; assert "PySide6" not in sys.modules' &&
        ../venv/bin/python -m microthermo validate --suite quick
      else
        ../venv/bin/python -m pip install --disable-pip-version-check \
          -r requirements-accel-lock.txt &&
        ../venv/bin/python -m pip install --disable-pip-version-check \
          '.[gui,accel]' &&
        QT_QPA_PLATFORM=offscreen ../venv/bin/python -m unittest \
          discover -s tests -q &&
        ../venv/bin/python -m microthermo validate --suite scientific
      fi
    ) &&
    "$batch/$label/venv/bin/python" -m pip freeze > "$batch/$label.packages.txt"
  } > "$batch/$label.stdout" 2> "$batch/$label.stderr"; then
    code=0
  else
    code=$?
  fi
  printf '%s\n' "$code" > "$batch/$label.exit"
  return "$code"
}

run_job headless & headless_pid=$!
run_job gui_accel & gui_pid=$!
headless_code=0
gui_code=0
wait "$headless_pid" || headless_code=$?
wait "$gui_pid" || gui_code=$?
printf 'headless exit: %s\ngui_accel exit: %s\n' "$headless_code" "$gui_code"
printf 'Give me this batch directory for review: %s\n' "$batch"
if (( headless_code != 0 || gui_code != 0 )); then exit 1; fi
