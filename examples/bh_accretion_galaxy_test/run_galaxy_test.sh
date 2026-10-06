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

rm -rf outputs/galaxy_env_shell
mkdir -p outputs/galaxy_env_shell logs

"$MPIEXEC" -np "$NTASKS" "$EXEC" param_galaxy.txt </dev/null 2>&1 | tee logs/galaxy_env_shell.log

if [[ ! -f outputs/galaxy_env_shell/end ]]; then
  echo "ERROR: galaxy run did not create outputs/galaxy_env_shell/end" >&2
  exit 1
fi

python3 analyze_galaxy_accretion.py
