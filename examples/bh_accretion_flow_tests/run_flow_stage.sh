#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
STAGE="${1:-diagnostic}"
EXEC="${EXEC:-$ROOT/ArepoBHFlows}"
MPIEXEC="${MPIEXEC:-mpirun}"
KEEP_RESTARTS="${KEEP_RESTARTS:-0}"

case "$STAGE" in
  diagnostic|evolved) ;;
  *) echo "Usage: $0 diagnostic|evolved" >&2; exit 2 ;;
esac

if [[ ! -x "$EXEC" ]]; then
  echo "Missing executable: $EXEC" >&2
  echo "Run ./build_flow_tests.sh first." >&2
  exit 1
fi

cd "$HERE"
python3 generate_run_files.py
python3 verify_ics.py
python3 verify_run_files.py

python3 - "$STAGE" <<'PY' |
import json, sys
stage=sys.argv[1]
m=json.load(open("run_manifest.json"))
for r in m["runs"]:
    if r["stage"] != stage:
        continue
    print("\t".join([
        r["parameter_file"], r["output_dir"], r["log_file"],
        str(r["recommended_ntasks"]), r["family"], r["tag"], r["level"]
    ]))
PY
while IFS=$'\t' read -r param out log ntasks family tag level
do
  echo
  echo "================================================================"
  echo " stage=$STAGE family=$family case=$tag level=$level ntasks=$ntasks"
  echo " param=$param"
  echo "================================================================"

  rm -rf "$out"
  mkdir -p "$out" "$(dirname "$log")"
  rm -f "$log"

  "$MPIEXEC" -np "$ntasks" "$EXEC" "$param" 2>&1 | tee "$log"

  if [[ ! -f "$out/end" ]]; then
    echo "ERROR: run did not create $out/end" >&2
    exit 1
  fi
  if ! grep -q 'BH_BENCHMARK_ALL:' "$log"; then
    echo "ERROR: missing BH_BENCHMARK_ALL diagnostic in $log" >&2
    exit 1
  fi
  if grep -Eqi 'terminate|segmentation fault|mpi_abort|floating point exception|(^|[^[:alpha:]])nan([^[:alpha:]]|$)' "$log"; then
    echo "ERROR: failure marker found in $log" >&2
    exit 1
  fi

  if [[ "$KEEP_RESTARTS" == "0" ]]; then
    find "$out" -maxdepth 1 -type f -name 'restart*' -delete
  fi

  grep -m1 'BH_BENCHMARK_ALL:' "$log"
done
