#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
CHECKPOINT="${CHECKPOINT:-$HERE/output/snap_001.hdf5}"
PARAM="$HERE/param_jet.txt"
OUT="$HERE/output_jet"
NTASKS="${NTASKS:-16}"

if [[ ! -x "$ROOT/ArepoSeedTest" ]]; then
  echo "Missing $ROOT/ArepoSeedTest; run ./build.sh first." >&2
  exit 1
fi

if [[ ! -f "$CHECKPOINT" ]]; then
  echo "Missing pre-seed checkpoint: $CHECKPOINT" >&2
  echo "Set CHECKPOINT=/path/to/preseed snapshot if needed." >&2
  exit 1
fi

python3 - "$CHECKPOINT" <<'PY'
import h5py
import sys
p=sys.argv[1]
with h5py.File(p,"r") as f:
    a=float(f["Header"].attrs["Time"])
    z=float(f["Header"].attrs["Redshift"])
    nbh=len(f["PartType5"]["ParticleIDs"]) if "PartType5" in f else 0
if nbh != 0:
    raise SystemExit(f"Checkpoint {p} already contains {nbh} BH(s); use the BH-free pre-seed snapshot.")
if not z > 47.7459:
    raise SystemExit(f"Checkpoint z={z} is not earlier than the deterministic seed event.")
print(f"Using BH-free checkpoint: a={a:.9f}, z={z:.6f}, N_BH={nbh}")
PY

cp "$HERE/param.txt" "$PARAM"
sed -i 's|OutputDir                               ./output|OutputDir                               ./output_jet|' "$PARAM"
sed -i 's|BHFreeFallA                             0.0|BHFreeFallA                             1.0e-3|' "$PARAM"
# This helper is a regression for an earlier iteration; disable Iteration-9
# internal accuracy limiting so its historical runtime/sampling stays stable.
sed -i 's|BHInternalTimestepFactor                1.0|BHInternalTimestepFactor                1.0e30|' "$PARAM"

# Jet-only Iteration-8 regression.  The large feedback aperture gives the
# narrow 15-degree bicone enough cells in this intentionally coarse smoke box.
# Wind release is suppressed so changes to gas kinetic energy are attributable
# to the jet channel alone.  Cmin=0.10 exercises the analytic jet-axis memory.
sed -i 's|BHFeedbackRadius                        0.20|BHFeedbackRadius                        1.00|' "$PARAM"
sed -i 's|BHWindBurstFactor                       0.01|BHWindBurstFactor                       1.0e30|' "$PARAM"
sed -i 's|BHJetBurstFactor                        0.10|BHJetBurstFactor                        0.10|' "$PARAM"
sed -i 's|BHMinActiveTargetMassFrac               0.5|BHMinActiveTargetMassFrac               0.01|' "$PARAM"
sed -i 's|BHMinCoherence                          0.30|BHMinCoherence                          0.10|' "$PARAM"

rm -rf "$OUT"
mkdir -p "$OUT"
cp "$CHECKPOINT" "$OUT/snap_001.hdf5"
rm -f "$HERE/run_jet.log"

cd "$HERE"
mpirun -np "$NTASKS" "$ROOT/ArepoSeedTest" "$PARAM" 2 1 2>&1 | tee "$HERE/run_jet.log"
