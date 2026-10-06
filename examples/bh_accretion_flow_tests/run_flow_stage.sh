#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
STAGE="${1:-diagnostic}"
EXEC="${EXEC:-$ROOT/ArepoBHFlows}"
MPIEXEC="${MPIEXEC:-mpirun}"
KEEP_RESTARTS="${KEEP_RESTARTS:-0}"
FORCE_RERUN="${FORCE_RERUN:-0}"

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

mapfile -t RUN_RECORDS < <(
python3 - "$STAGE" <<'PY'
import json
import sys

stage = sys.argv[1]
with open("run_manifest.json", "r", encoding="utf-8") as f:
    manifest = json.load(f)

for run in manifest["runs"]:
    if run["stage"] != stage:
        continue
    print("\t".join([
        run["parameter_file"],
        run["output_dir"],
        run["log_file"],
        str(run["recommended_ntasks"]),
        run["family"],
        run["tag"],
        run["level"],
    ]))
PY
)

TOTAL="${#RUN_RECORDS[@]}"
if [[ "$TOTAL" -eq 0 ]]; then
  echo "No runs found for stage '$STAGE'." >&2
  exit 1
fi

echo
echo "Prepared $TOTAL runs for stage=$STAGE"

INDEX=0
for record in "${RUN_RECORDS[@]}"
do
  INDEX=$((INDEX + 1))

  param=""
  out=""
  log=""
  ntasks=""
  family=""
  tag=""
  level=""
  IFS=$'\t' read -r param out log ntasks family tag level <<< "$record"

  if [[ -z "$param" || -z "$out" || -z "$log" || -z "$ntasks" ||
        -z "$family" || -z "$tag" || -z "$level" ]]; then
    echo "ERROR: malformed run record: $record" >&2
    exit 1
  fi

  echo
  echo "================================================================"
  echo " run=$INDEX/$TOTAL stage=$STAGE family=$family case=$tag level=$level ntasks=$ntasks"
  echo " param=$param"
  echo "================================================================"

  if [[ "$FORCE_RERUN" == "0" && -f "$out/end" && -f "$log" ]] &&
     grep -q 'BH_BENCHMARK_ALL:' "$log" &&
     ! grep -Eqi 'terminate|segmentation fault|mpi_abort|floating point exception|(^|[^[:alpha:]])nan([^[:alpha:]]|$)' "$log"
  then
    echo " SKIP: already completed successfully"
    grep -m1 'BH_BENCHMARK_ALL:' "$log"
    continue
  fi

  rm -rf "$out"
  mkdir -p "$out" "$(dirname "$log")"
  rm -f "$log"

  # The run list is already stored in RUN_RECORDS, and MPI receives /dev/null
  # on stdin so it cannot consume runner control data.
  "$MPIEXEC" -np "$ntasks" "$EXEC" "$param" </dev/null 2>&1 | tee "$log"

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

echo
echo "Completed stage=$STAGE: $TOTAL/$TOTAL manifest entries processed."
