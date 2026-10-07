#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

EXEC="${EXEC:-$ROOT/ArepoRarePeak}"
MPIEXEC="${MPIEXEC:-mpirun}"
NTASKS="${NTASKS:-4}"
COMMON_RESTART_DIR="${COMMON_RESTART_DIR:-$HERE/restart_archive/common_seed_spare/restartfiles}"
PARAM_DIR="$HERE/params_native"
LOG_DIR="$HERE/logs"
BRANCHES="${BRANCHES:-}"

available_branches()
{
  for p in "$PARAM_DIR"/*.txt; do
    [[ -e "$p" ]] || continue
    basename "$p" .txt
  done
}

if [[ ! -x "$EXEC" ]]; then
  echo "Missing executable: $EXEC. Run ./build.sh first." >&2
  exit 1
fi

if [[ ! -d "$COMMON_RESTART_DIR" ]]; then
  echo "Missing common native restart directory: $COMMON_RESTART_DIR" >&2
  echo "Expected the untouched post-seed z=20.8837 restart checkpoint." >&2
  exit 1
fi

if [[ -z "$BRANCHES" ]]; then
  echo "No branches selected. Set BRANCHES explicitly; nothing was run."
  echo
  echo "Available branches:"
  available_branches | sed 's/^/  /'
  echo
  echo 'Example:'
  echo '  BRANCHES="convj_shell_ffr__macer" NTASKS=4 ./run_native_restart_suite.sh'
  exit 0
fi

mkdir -p "$LOG_DIR"

for name in $BRANCHES; do
  param="$PARAM_DIR/$name.txt"
  if [[ ! -s "$param" ]]; then
    echo "Unknown branch '$name'. Available:" >&2
    available_branches | sed 's/^/  /' >&2
    exit 1
  fi

  out_rel="$(awk '$1=="OutputDir" {print $2; exit}' "$param")"
  if [[ -z "$out_rel" ]]; then
    echo "No OutputDir in $param" >&2
    exit 1
  fi
  if [[ "$out_rel" == ./* ]]; then
    out="$HERE/${out_rel#./}"
  else
    out="$out_rel"
  fi

  log="$LOG_DIR/$name.log"

  # Science branches must always begin from a pristine copy of the common
  # post-seed native restart. Never delete or overwrite an existing result.
  if [[ -e "$out/restartfiles" || -e "$out/end" || -s "$log" ]]; then
    echo "Refusing to overwrite existing branch state for '$name'." >&2
    echo "  output: $out" >&2
    echo "  log:    $log" >&2
    echo "Choose a new branch/output name or archive the existing result manually." >&2
    exit 1
  fi

  mkdir -p "$out"
  cp -a --reflink=auto "$COMMON_RESTART_DIR" "$out/restartfiles"

  echo
  echo "================================================================"
  echo " native-restart branch: $name"
  echo " checkpoint: z=20.8837, BH ID=1000120011, M_BH=1e5 Msun"
  echo " MPI ranks: $NTASKS  (must match the native checkpoint decomposition)"
  echo " output: $out"
  echo "================================================================"

  (
    cd "$HERE"
    "$MPIEXEC" -np "$NTASKS" "$EXEC" "$param" 1 </dev/null 2>&1 | tee "$log"
  )

  if [[ ! -f "$out/end" ]]; then
    echo "ERROR: branch '$name' did not finish cleanly." >&2
    exit 1
  fi
done
