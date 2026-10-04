#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
SOURCE="${SOURCE:-$HERE/output_dm/snap_003.hdf5}"
PARAM="$HERE/param_restart_guard.txt"
OUT="$HERE/output_restart_guard"
LOG="$HERE/run_restart_guard.log"
NTASKS="${NTASKS:-4}"

if [[ ! -x "$ROOT/ArepoSeedTest" ]]; then
  echo "Missing $ROOT/ArepoSeedTest; building the smoke-test executable now." >&2
  "$HERE/build.sh"
fi
if [[ ! -x "$ROOT/ArepoSeedTest" ]]; then
  echo "Build completed without an executable at $ROOT/ArepoSeedTest." >&2
  exit 1
fi
if [[ ! -f "$SOURCE" ]]; then
  echo "Missing evolved one-BH snapshot: $SOURCE" >&2
  exit 1
fi

cp "$HERE/param.txt" "$PARAM"
sed -i 's|OutputDir                               ./output|OutputDir                               ./output_restart_guard|' "$PARAM"

rm -rf "$OUT"
mkdir -p "$OUT"
cp "$SOURCE" "$OUT/snap_003.hdf5"
rm -f "$LOG"

cd "$HERE"
set +e
mpirun -np "$NTASKS" "$ROOT/ArepoSeedTest" "$PARAM" 2 3 >"$LOG" 2>&1
status=$?
set -e

if [[ $status -eq 0 ]]; then
  echo "FAIL: evolved Type-5 RestartFlag=2 unexpectedly succeeded" >&2
  exit 1
fi

if ! grep -q "RestartFlag=2 from a snapshot containing" "$LOG"; then
  echo "FAIL: RestartFlag=2 failed, but not through the explicit FFR state-safety guard" >&2
  tail -n 80 "$LOG" >&2
  exit 1
fi

echo "PASS: evolved Type-5 RestartFlag=2 is explicitly rejected"
echo "      native RestartFlag=1 remains the supported evolved FFR restart path"
