#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

EXEC="${EXEC:-$ROOT/ArepoRarePeak}"
MPIEXEC="${MPIEXEC:-mpirun}"
NTASKS="${NTASKS:-4}"
MIN_FREE_GIB="${MIN_FREE_GIB:-20}"
COMMON_RESTART_DIR="${COMMON_RESTART_DIR:-$HERE/restart_archive/common_seed_spare/restartfiles}"
PARAM_DIR="$HERE/params_native"
LOG_DIR="$HERE/logs"
BRANCHES="${BRANCHES:-}"
RESUME="${RESUME:-0}"

available_branches()
{
  local p
  for p in "$PARAM_DIR"/*.txt; do
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

restart_files=("$COMMON_RESTART_DIR"/restart.[0-9]*)
restart_tasks=${#restart_files[@]}

if [[ "$restart_tasks" -lt 1 ]]; then
  echo "No native restart.N files found in $COMMON_RESTART_DIR" >&2
  exit 1
fi

if [[ "$NTASKS" -ne "$restart_tasks" ]]; then
  echo "Native restart rank mismatch: checkpoint has $restart_tasks rank-local files, NTASKS=$NTASKS." >&2
  echo "AREPO native RestartFlag=1 requires the same MPI rank count as the writer." >&2
  exit 1
fi

avail_kib="$(df -Pk "$HERE" | awk 'NR==2 {print $4}')"
min_kib=$((MIN_FREE_GIB * 1024 * 1024))

if [[ "$avail_kib" -lt "$min_kib" ]]; then
  echo "Insufficient free space for a science branch: need >=${MIN_FREE_GIB} GiB free." >&2
  df -h "$HERE" >&2
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

  if [[ "$RESUME" == "1" ]]; then
    if [[ ! -d "$out/restartfiles" ]]; then
      echo "Cannot resume '$name': missing $out/restartfiles" >&2
      exit 1
    fi
    if [[ -e "$out/end" ]]; then
      echo "Cannot resume '$name': branch already has an end file." >&2
      exit 1
    fi
    # A deliberate AREPO stop leaves this sentinel behind. Remove only this
    # branch-local stop request so RestartFlag=1 can advance again.
    rm -f -- "$out/stop"
  else
    # Every fresh science branch starts from an independent copy of the
    # pristine post-seed native checkpoint. Existing branch data are never
    # removed or overwritten.
    if [[ -e "$out/restartfiles" || -e "$out/end" || -s "$log" ]]; then
      echo "Refusing to overwrite existing branch state for '$name'." >&2
      echo "  output: $out" >&2
      echo "  log:    $log" >&2
      echo "Use RESUME=1 for a verified native restart, or archive the existing result manually." >&2
      exit 1
    fi

    mkdir -p "$out"
    cp -a --reflink=auto "$COMMON_RESTART_DIR" "$out/restartfiles"
  fi

  echo
  echo "================================================================"
  echo " native-restart branch: $name"
  if [[ "$RESUME" == "1" ]]; then
    echo " mode: resume existing verified native restart"
  else
    echo " checkpoint: z=20.8837, BH ID=1000120011, M_BH=1e5 Msun"
  fi
  echo " MPI ranks: $NTASKS"
  echo " free-space floor: ${MIN_FREE_GIB} GiB"
  echo " output: $out"
  echo "================================================================"

  if [[ "$RESUME" == "1" && -f "$log" ]]; then
    run_log_start=$(( $(wc -l < "$log") + 1 ))
  else
    run_log_start=1
  fi

  (
    cd "$HERE"
    if [[ "$RESUME" == "1" ]]; then
      "$MPIEXEC" -np "$NTASKS" "$EXEC" "$param" 1 </dev/null 2>&1 | tee -a "$log"
    else
      "$MPIEXEC" -np "$NTASKS" "$EXEC" "$param" 1 </dev/null 2>&1 | tee "$log"
    fi
  )

  if [[ -f "$out/end" ]]; then
    echo "Branch '$name' reached its configured final time."
  elif [[ -d "$out/restartfiles" ]] &&
       tail -n +"$run_log_start" "$log" | grep -q "stop-file detected. stopping." &&
       tail -n +"$run_log_start" "$log" | grep -q "All restart files written successfully"; then
    echo "Branch '$name' stopped cleanly with verified native restart files."
  else
    echo "ERROR: branch '$name' exited without an end file or a verified deliberate-stop restart." >&2
    exit 1
  fi
done
