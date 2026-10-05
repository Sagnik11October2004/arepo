#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

SOURCE="${SOURCE:-$HERE/accretion_source.hdf5}"
DA="${DA:-8.0e-8}"
STEP_FRAC="${STEP_FRAC:-0.45}"
NTASKS="${NTASKS:-16}"
EXEC="${EXEC:-$ROOT/ArepoFFRSnapshotTest}"
ESTIMATOR_EDD_FACTOR="${ESTIMATOR_EDD_FACTOR:-1.0e30}"
RSINK_TOL="${RSINK_TOL:-1.0e-3}"
RAW_STABILITY_TOL="${RAW_STABILITY_TOL:-5.0e-2}"

if [[ ! -x "$EXEC" ]]; then
  echo "Missing $EXEC; run ./build_snapshot_ic.sh first." >&2
  exit 1
fi
if [[ ! -f "$SOURCE" ]]; then
  echo "Missing durable source snapshot: $SOURCE" >&2
  echo "Create it with: ./prepare_accretion_source.sh /path/to/source_snapshot.hdf5" >&2
  exit 1
fi
if [[ -L "$SOURCE" ]]; then
  echo "Refusing symlinked accretion source $SOURCE; use prepare_accretion_source.sh to make a standalone copy." >&2
  exit 1
fi
case "$SOURCE" in
  *.hdf5) ;;
  *) echo "Accretion source must end in .hdf5: $SOURCE" >&2; exit 1 ;;
esac

SOURCE_ABS="$(python3 - "$SOURCE" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).resolve())
PY
)"
SOURCE_STEM="${SOURCE_ABS%.hdf5}"

python3 - "$SOURCE_ABS" <<'PY'
import h5py, sys
p=sys.argv[1]
with h5py.File(p,"r") as f:
    a=float(f["Header"].attrs["Time"])
    z=float(f["Header"].attrs["Redshift"])
    nbh=len(f["PartType5"]["ParticleIDs"]) if "PartType5" in f else 0
    if nbh != 1:
        raise SystemExit(f"Expected exactly one BH in {p}; found {nbh}")
    disk=f["PartType5"].get("BH_DiskMass")
    wind=f["PartType5"].get("BH_WindBufferMass")
    if disk is not None and float(disk[0]) != 0.0:
        raise SystemExit(f"Direct-backend smoke test requires zero BH_DiskMass; found {float(disk[0])}")
    if wind is not None and float(wind[0]) != 0.0:
        raise SystemExit(f"Direct-backend smoke test requires zero BH_WindBufferMass; found {float(wind[0])}")
print(f"Using durable source: {p} a={a:.15f} z={z:.9f} N_BH={nbh}")
PY

read A0 AMAX DLOGA MAXSTEP <<<"$(python3 - "$SOURCE_ABS" "$DA" "$STEP_FRAC" <<'PY'
import h5py, math, sys
p=sys.argv[1]
da=float(sys.argv[2])
step_frac=float(sys.argv[3])
if not da > 0:
    raise SystemExit(f"DA must be positive, got {da}")
if not (0.0 < step_frac < 1.0):
    raise SystemExit(f"STEP_FRAC must lie in (0,1), got {step_frac}")
with h5py.File(p,"r") as f:
    a=float(f["Header"].attrs["Time"])
amax=a+da
dloga=math.log(amax/a)
maxstep=step_frac*dloga
if not (0.0 < maxstep < dloga):
    raise SystemExit(f"invalid comoving MaxSizeTimestep={maxstep} for dloga={dloga}")
print(f"{a:.15f} {amax:.15f} {dloga:.15e} {maxstep:.15e}")
PY
)"

echo "TimeBegin       = $A0"
echo "TimeMax         = $AMAX"
echo "Delta ln(a)     = $DLOGA"
echo "STEP_FRAC       = $STEP_FRAC"
echo "MaxSizeTimestep = $MAXSTEP"
echo "Estimator-only Eddington factor = $ESTIMATOR_EDD_FACTOR"

run_one()
{
  local name="$1"
  local model="$2"
  local param="$HERE/param_acc_${name}.txt"
  local out="$HERE/output_acc_${name}"
  local log="$HERE/run_acc_${name}.log"

  cp "$HERE/param.txt" "$param"

  python3 - "$param" "$A0" "$AMAX" "$MAXSTEP" "$model" "$out" "$SOURCE_STEM" "$ESTIMATOR_EDD_FACTOR" <<'PY'
from pathlib import Path
import sys

p=Path(sys.argv[1])
a0,amax,maxstep=sys.argv[2:5]
model=sys.argv[5]
out=sys.argv[6]
source_stem=sys.argv[7]
edd_factor=sys.argv[8]

updates={
    "InitCondFile":source_stem,
    "OutputDir":out,
    "TimeBegin":a0,
    "TimeMax":amax,
    "MaxSizeTimestep":maxstep,
    "OutputListOn":"0",
    # In comoving mode AREPO explicitly requires TimeBetSnapshot > 1.
    # First snapshot time is beyond this short run, so no science snapshot is
    # requested by the smoke harness.
    "TimeBetSnapshot":"2.0",
    "TimeOfFirstSnapshot":"1.0",
    "BHSeedMinRedshift":"100.0",
    "BHFreeFallA":"1.0e-3",
    "BHInternalTimestepFactor":"1.0e30",
    "BHDMNeighbours":"2",
    "BHBenchmarkAccretionModel":model,
    "BHBenchmarkAccretionTarget":"1",
    "BHBenchmarkFeedbackModel":"0",
    # Estimator-only phase: prevent the TNG backend cap from hiding differences
    # in the raw resolved prescriptions. Fiducial TNG runs restore this to 1.
    "BHBenchmarkEddingtonFactor":edd_factor,
}

lines=p.read_text().splitlines()
seen=set()
out_lines=[]
for line in lines:
    stripped=line.strip()
    if not stripped or stripped.startswith("%"):
        out_lines.append(line)
        continue
    key=stripped.split()[0]
    if key in updates:
        out_lines.append(f"{key:<40s}{updates[key]}")
        seen.add(key)
    else:
        out_lines.append(line)

missing=set(updates)-seen
if missing:
    raise SystemExit("Missing parameter keys: "+", ".join(sorted(missing)))

p.write_text("\n".join(out_lines)+"\n")
PY

  rm -rf "$out"
  mkdir -p "$out"
  rm -f "$log"

  echo
  echo "=================================================="
  echo " Running $name   model=$model direct + NONE"
  echo "=================================================="

  (
    cd "$HERE"
    mpirun -np "$NTASKS" "$EXEC" "$(basename "$param")" 2>&1 | tee "$log"
  )
}

run_one tng_bondi 0
run_one boosted_bondi 1
run_one am_bondi 2
run_one ffr 3

python3 "$HERE/verify_accretion_modes.py" \
  --rsink-tol "$RSINK_TOL" \
  --raw-stability-tol "$RAW_STABILITY_TOL" \
  --expected-edd-factor "$ESTIMATOR_EDD_FACTOR"

echo
echo "================ first diagnostic from each model ================"
for name in tng_bondi boosted_bondi am_bondi ffr
do
  printf "%-18s " "$name"
  grep -m1 "BH_BENCHMARK: accretion" "$HERE/run_acc_${name}.log"
done
