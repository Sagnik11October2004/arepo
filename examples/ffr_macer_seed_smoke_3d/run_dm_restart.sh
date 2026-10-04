#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
CHECKPOINT="${CHECKPOINT:-$HERE/output/snap_001.hdf5}"
PARAM="$HERE/param_dm.txt"
OUT="$HERE/output_dm"
NTASKS="${NTASKS:-16}"

if [[ ! -x "$ROOT/ArepoSeedTest" ]]; then
  echo "Missing $ROOT/ArepoSeedTest; building the smoke-test executable now." >&2
  "$HERE/build.sh"
fi
if [[ ! -x "$ROOT/ArepoSeedTest" ]]; then
  echo "Build completed without an executable at $ROOT/ArepoSeedTest." >&2
  exit 1
fi

if [[ ! -f "$CHECKPOINT" ]]; then
  echo "Missing pre-seed checkpoint: $CHECKPOINT" >&2
  exit 1
fi

python3 - "$CHECKPOINT" <<'PY'
import h5py, sys
p=sys.argv[1]
with h5py.File(p,"r") as f:
    a=float(f["Header"].attrs["Time"])
    z=float(f["Header"].attrs["Redshift"])
    nbh=len(f["PartType5"]["ParticleIDs"]) if "PartType5" in f else 0
if nbh != 0:
    raise SystemExit(f"Checkpoint {p} already contains {nbh} BH(s)")
if not z > 47.7459:
    raise SystemExit(f"Checkpoint z={z} is not earlier than the deterministic seed event")
print(f"Using BH-free checkpoint: a={a:.9f}, z={z:.6f}, N_BH={nbh}")
PY

cp "$HERE/param.txt" "$PARAM"
sed -i 's|OutputDir                               ./output|OutputDir                               ./output_dm|' "$PARAM"

# Iteration-10 isolation:
# - canonical BHFreeFallA=0 means no reservoir/feedback energy is generated;
# - central binding is disabled, so positive Eth values can only come from
#   the cached local DM velocity dispersion;
# - internal timestep limiting is disabled for historical smoke cadence.
sed -i 's|BHUseCentralBindingTerm                 1|BHUseCentralBindingTerm                 0|' "$PARAM"
sed -i 's|BHInternalTimestepFactor                1.0|BHInternalTimestepFactor                1.0e30|' "$PARAM"

rm -rf "$OUT"
mkdir -p "$OUT"
cp "$CHECKPOINT" "$OUT/snap_001.hdf5"
rm -f "$HERE/run_dm.log"

cd "$HERE"
mpirun -np "$NTASKS" "$ROOT/ArepoSeedTest" "$PARAM" 2 1 2>&1 | tee "$HERE/run_dm.log"
