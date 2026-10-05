#!/usr/bin/env bash
# TEMPORARY launcher: delete after the archive and replay have been reviewed.
# The archive and logs under runs/carnot_triangles_500_seed123_3cycles are retained.
set -u
cd "$(dirname "$0")/.." || exit 1

out=runs/carnot_triangles_500_seed123_3cycles
archive="$out/archive"

case "${1:-}" in
  record)
    if [[ -e "$archive" ]]; then
      echo "Archive already exists: $archive" >&2
      echo 'Choose a new output path or remove the incomplete archive before retrying.' >&2
      exit 2
    fi
    mkdir -p "$out"
    if .venv/bin/python -m microthermo precalculate \
        --preset carnot_triangles --particles 500 --seed 123 \
        --shaft-speed 0.15 --duration 125.66370614359172 \
        --fps 60 --chunk-frames 256 --output "$archive" \
        >"$out/precalculate.stdout" 2>"$out/precalculate.stderr"; then
      code=0
    else
      code=$?
    fi
    printf '%s\n' "$code" > "$out/precalculate.exit"
    printf 'Precalculation exit: %s; output: %s\n' "$code" "$out"
    exit "$code"
    ;;
  replay)
    if [[ ! -f "$archive/manifest.json" || ! -f "$out/precalculate.exit" ]]; then
      echo "No complete archive at $archive; finish the record step first." >&2
      exit 2
    fi
    if [[ $(<"$out/precalculate.exit") != 0 ]]; then
      echo "No complete archive at $archive; finish the record step first." >&2
      exit 2
    fi
    echo 'Opening replay window; press Play to start.'
    if .venv/bin/python -m microthermo replay "$archive" \
        >"$out/replay.stdout" 2>"$out/replay.stderr"; then
      code=0
    else
      code=$?
    fi
    printf '%s\n' "$code" > "$out/replay.exit"
    exit "$code"
    ;;
  *)
    echo "Usage: bash $0 {record|replay}" >&2
    exit 2
    ;;
esac
