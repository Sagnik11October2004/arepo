#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
NTASKS="${NTASKS:-4}"

if [[ ! -x "$ROOT/ArepoSeedTest" ]]; then
  echo "Missing $ROOT/ArepoSeedTest. Run ./build.sh first." >&2
  exit 1
fi
if [[ ! -s "$HERE/ics.hdf5" ]]; then
  echo "Missing $HERE/ics.hdf5. Run ./generate_ics.sh first." >&2
  exit 1
fi

cd "$HERE"
rm -rf output
mkdir -p output
rm -f run.log

mpirun -np "$NTASKS" "$ROOT/ArepoSeedTest" param.txt 2>&1 | tee run.log
