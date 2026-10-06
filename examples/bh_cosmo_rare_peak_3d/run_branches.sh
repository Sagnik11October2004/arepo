#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
EXEC="${EXEC:-$ROOT/ArepoRarePeak}"
MPIEXEC="${MPIEXEC:-mpirun}"
NTASKS="${NTASKS:-8}"
MODELS="${MODELS:-0 1 2 3 4 5}"
FEEDBACKS="${FEEDBACKS:-none tng macer}"
ZEND="${ZEND:-15}"
NONE_TARGET="${NONE_TARGET:-direct}"

if [[ ! -x "$EXEC" ]]; then
  echo "Missing executable: $EXEC. Run ./build.sh first." >&2
  exit 1
fi
if [[ ! -s "$HERE/ics/common_z20_seeded.hdf5" ]]; then
  echo "Missing common z=20 seeded IC. Run ./run_to_z20.sh first." >&2
  exit 1
fi

cd "$HERE"
python3 prepare_branch_params.py --z-end "$ZEND" --none-target "$NONE_TARGET"
mkdir -p logs outputs

model_name()
{
  case "$1" in
    0) echo tng_bondi ;;
    1) echo boosted_bondi ;;
    2) echo am_bondi ;;
    3) echo ffr_volume ;;
    4) echo ffr_shell ;;
    5) echo convj_shell_ffr ;;
    *) echo "invalid model id $1" >&2; return 1 ;;
  esac
}

for model in $MODELS; do
  name="$(model_name "$model")"
  for fb in $FEEDBACKS; do
    case "$fb" in
      none|tng|macer) ;;
      *) echo "invalid feedback '$fb' (use none, tng, macer)" >&2; exit 1 ;;
    esac

    param="params/${name}__${fb}.txt"
    out="outputs/${name}__${fb}"
    log="logs/${name}__${fb}.log"

    rm -rf "$out"
    mkdir -p "$out"
    rm -f "$log"

    echo
    echo "================================================================"
    echo " model=$model ($name)   feedback=$fb   z=20 -> $ZEND"
    echo "================================================================"
    "$MPIEXEC" -np "$NTASKS" "$EXEC" "$param" </dev/null 2>&1 | tee "$log"

    if [[ ! -f "$out/end" ]]; then
      echo "ERROR: branch $name/$fb did not finish cleanly." >&2
      exit 1
    fi
  done
done
