#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEP_FRACS="${STEP_FRACS:-0.45 0.225 0.1125 0.05625}"
OUTROOT="${OUTROOT:-$HERE/sink_timestep_convergence}"
NTASKS="${NTASKS:-16}"
DA="${DA:-8.0e-8}"

rm -rf "$OUTROOT"
mkdir -p "$OUTROOT"

for frac in $STEP_FRACS
do
  tag="$(python3 - "$frac" <<'PY'
import sys
print(("%g" % float(sys.argv[1])).replace(".","p"))
PY
)"
  echo
  echo "################ sink convergence STEP_FRAC=$frac ################"
  STEP_FRAC="$frac" NTASKS="$NTASKS" DA="$DA" RSINK_TOL=5.0e-2 \
    "$HERE/run_accretion_modes.sh"

  dst="$OUTROOT/step_$tag"
  mkdir -p "$dst"
  cp "$HERE"/run_acc_*.log "$dst"/
  cp "$HERE"/param_acc_*.txt "$dst"/

  # These short fresh-IC runs can leave large native restart files even though
  # no science snapshot is requested. The archived logs/params contain the
  # convergence diagnostics; remove run products between timestep levels.
  rm -rf "$HERE"/output_acc_tng_bondi "$HERE"/output_acc_boosted_bondi \
         "$HERE"/output_acc_am_bondi "$HERE"/output_acc_ffr
done

python3 "$HERE/verify_sink_timestep_convergence.py" \
  --root "$OUTROOT" --step-fracs $STEP_FRACS
