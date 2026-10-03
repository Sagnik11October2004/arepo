#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MUSIC_BIN="${MUSIC_BIN:-$HOME/Research/MUSIC2/build/MUSIC}"

if [[ ! -x "$MUSIC_BIN" ]]; then
  echo "MUSIC executable not found: $MUSIC_BIN" >&2
  echo "Set MUSIC_BIN=/absolute/path/to/MUSIC and rerun." >&2
  exit 1
fi

cd "$HERE"
rm -f ics.hdf5
"$MUSIC_BIN" music.conf

if [[ ! -s ics.hdf5 ]]; then
  echo "MUSIC did not create $HERE/ics.hdf5" >&2
  exit 1
fi

echo "Created: $HERE/ics.hdf5"
