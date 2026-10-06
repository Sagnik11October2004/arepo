#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
EXEC="${EXEC:-$ROOT/ArepoBHFlows}"
MPIEXEC="${MPIEXEC:-mpirun}"
NTASKS="${NTASKS:-8}"

cd "$HERE"
python3 prepare_galaxy_test.py

if [[ ! -x "$EXEC" ]]; then
  echo "Missing executable: $EXEC" >&2
  echo "Build first with ../bh_accretion_flow_tests/build_flow_tests.sh" >&2
  exit 1
fi

rm -rf outputs logs
mkdir -p outputs logs

echo "=== Stage 1: common BH-off settling run (0-20 Myr) ==="
"$MPIEXEC" -np "$NTASKS" "$EXEC" param_settle.txt </dev/null 2>&1 | tee logs/settle.log
if [[ ! -f outputs/settle/end ]]; then
  echo "ERROR: settle run did not create outputs/settle/end" >&2
  exit 1
fi

SETTLED_SNAPSHOT="$(ls outputs/settle/snap_*.hdf5 2>/dev/null | sort | tail -n 1)"
if [[ -z "$SETTLED_SNAPSHOT" ]]; then
  echo "ERROR: no settled HDF5 snapshot found" >&2
  exit 1
fi
echo "Converting settled state into a fresh t=0 branch IC: $SETTLED_SNAPSHOT"
python3 make_settled_branch_ic.py "$SETTLED_SNAPSHOT"

MODELS=(
  "0 tng_bondi"
  "1 boosted_bondi"
  "2 am_bondi"
  "3 ffr_volume"
  "4 ffr_shell"
  "5 ffr_env"
)

for entry in "${MODELS[@]}"; do
  read -r model name <<< "$entry"
  param="params/model_${model}_${name}.txt"
  out="outputs/model_${model}_${name}"
  log="logs/model_${model}_${name}.log"
  rm -rf "$out"
  echo
  echo "=== Stage 2: model=$model ($name), common settled IC, 0-60 Myr branch ==="
  "$MPIEXEC" -np "$NTASKS" "$EXEC" "$param" </dev/null 2>&1 | tee "$log"
  if [[ ! -f "$out/end" ]]; then
    echo "ERROR: branch $model/$name did not create $out/end" >&2
    exit 1
  fi
done

python3 analyze_galaxy_accretion.py
