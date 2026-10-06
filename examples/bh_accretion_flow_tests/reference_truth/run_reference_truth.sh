#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
FLOW="$(cd "$HERE/.." && pwd)"
EXEC="\${EXEC:-$FLOW/../../ArepoBHFlows}"
MPIEXEC="\${MPIEXEC:-mpirun}"
FORCE_RERUN="\${FORCE_RERUN:-0}"

# Optional filters, e.g. FAMILY=bhl TAG=M5p0 LEVEL=HR6 ./run_reference_truth.sh
FAMILY_FILTER="\${FAMILY:-}"
TAG_FILTER="\${TAG:-}"
LEVEL_FILTER="\${LEVEL:-}"

cd "$FLOW"

python3 reference_truth/generate_reference_ics.py
python3 reference_truth/generate_reference_runs.py

mapfile -t RUN_RECORDS < <(
python3 - "$FAMILY_FILTER" "$TAG_FILTER" "$LEVEL_FILTER" <<'PY'
import json
import sys

family_filter, tag_filter, level_filter = sys.argv[1:]
m = json.load(open("reference_truth/reference_run_manifest.json"))

for r in m["runs"]:
    if family_filter and r["family"] != family_filter:
        continue
    if tag_filter and r["tag"] != tag_filter:
        continue
    if level_filter and r["level"] != level_filter:
        continue
    print("\t".join([
        r["parameter_file"],
        r["output_dir"],
        r["log_file"],
        str(r["recommended_ntasks"]),
        r["family"],
        r["tag"],
        r["level"],
    ]))
PY
)

TOTAL="\${#RUN_RECORDS[@]}"
if [[ "$TOTAL" -eq 0 ]]; then
  echo "No reference runs match the requested filters." >&2
  exit 1
fi

if [[ ! -x "$EXEC" ]]; then
  echo "Missing executable: $EXEC" >&2
  echo "Build the flow-test executable first with ./build_flow_tests.sh" >&2
  exit 1
fi

echo
echo "Prepared $TOTAL high-resolution reference runs."

INDEX=0
for record in "\${RUN_RECORDS[@]}"
do
  INDEX=$((INDEX + 1))
  IFS=$'\t' read -r param out log ntasks family tag level <<< "$record"

  echo
  echo "================================================================"
  echo " reference=$INDEX/$TOTAL family=$family case=$tag level=$level ntasks=$ntasks"
  echo " param=$param"
  echo "================================================================"

  if [[ "$FORCE_RERUN" == "0" && -f "$out/end" && -f "$log" ]] &&
     grep -q 'BH_BENCHMARK: mass-ledger' "$log" &&
     ! grep -Eqi 'terminate|segmentation fault|mpi_abort|floating point exception|(^|[^[:alpha:]])nan([^[:alpha:]]|$)' "$log"
  then
    echo " SKIP: already completed successfully"
    continue
  fi

  rm -rf "$out"
  mkdir -p "$out" "$(dirname "$log")"
  rm -f "$log"

  "$MPIEXEC" -np "$ntasks" "$EXEC" "$param" </dev/null 2>&1 | tee "$log"

  if [[ ! -f "$out/end" ]]; then
    echo "ERROR: reference run did not create $out/end" >&2
    exit 1
  fi

  if grep -Eqi 'terminate|segmentation fault|mpi_abort|floating point exception|(^|[^[:alpha:]])nan([^[:alpha:]]|$)' "$log"; then
    echo "ERROR: failure marker found in $log" >&2
    exit 1
  fi
done

echo
echo "Completed selected high-resolution reference runs."
