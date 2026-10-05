#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

SOURCE="${SOURCE:-$HERE/output_dm/snap_003.hdf5}"
DA="${DA:-8.0e-8}"
NTASKS="${NTASKS:-16}"
EXEC="${EXEC:-$ROOT/ArepoFFRSnapshotTest}"

if [[ ! -x "$EXEC" ]]; then
  echo "Missing $EXEC; run ./build_snapshot_ic.sh first." >&2
  exit 1
fi
if [[ ! -f "$SOURCE" ]]; then
  echo "Missing source snapshot: $SOURCE" >&2
  exit 1
fi

python3 - "$SOURCE" <<'PY'
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
print(f"Using source: a={a:.15f} z={z:.9f} N_BH={nbh}")
PY

ln -sfn "$SOURCE" "$HERE/accretion_source.hdf5"

read A0 AMAX MAXSTEP <<<"$(python3 - "$SOURCE" "$DA" <<'PY'
import h5py, math, sys
p=sys.argv[1]
da=float(sys.argv[2])
with h5py.File(p,"r") as f:
    a=float(f["Header"].attrs["Time"])
amax=a+da
dloga=math.log(amax/a)
maxstep=0.45*min(da,dloga)
print(f"{a:.15f} {amax:.15f} {maxstep:.15e}")
PY
)"

echo "TimeBegin       = $A0"
echo "TimeMax         = $AMAX"
echo "MaxSizeTimestep = $MAXSTEP"

run_one()
{
  local name="$1"
  local model="$2"
  local param="$HERE/param_acc_${name}.txt"
  local out="$HERE/output_acc_${name}"
  local log="$HERE/run_acc_${name}.log"

  cp "$HERE/param.txt" "$param"

  python3 - "$param" "$A0" "$AMAX" "$MAXSTEP" "$model" "$out" <<'PY'
from pathlib import Path
import sys

p=Path(sys.argv[1])
a0,amax,maxstep=sys.argv[2:5]
model=sys.argv[5]
out=sys.argv[6]

updates={
    "InitCondFile":"./accretion_source",
    "OutputDir":out,
    "TimeBegin":a0,
    "TimeMax":amax,
    "MaxSizeTimestep":maxstep,
    "OutputListOn":"0",
    # In comoving mode AREPO explicitly requires TimeBetSnapshot > 1.
    # First snapshot time is far beyond this very short smoke run, so no
    # snapshot is actually written.
    "TimeBetSnapshot":"2.0",
    "TimeOfFirstSnapshot":"1.0",
    "BHSeedMinRedshift":"100.0",
    "BHFreeFallA":"1.0e-3",
    "BHInternalTimestepFactor":"1.0e30",
    "BHDMNeighbours":"2",
    "BHBenchmarkAccretionModel":model,
    "BHBenchmarkAccretionTarget":"1",
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
  echo " Running $name   model=$model"
  echo "=================================================="

  (
    cd "$HERE"
    mpirun -np "$NTASKS" "$EXEC" "$(basename "$param")" 2>&1 | tee "$log"
  )

  echo
  echo "Accretion diagnostics for $name:"
  grep "BH_BENCHMARK: accretion" "$log" || {
    echo "ERROR: no BH_BENCHMARK accretion diagnostic was produced." >&2
    return 1
  }
}

run_one tng_bondi 0
run_one boosted_bondi 1
run_one am_bondi 2
run_one ffr 3

echo
echo "================ first diagnostic from each model ================"
for name in tng_bondi boosted_bondi am_bondi ffr
do
  printf "%-18s " "$name"
  grep -m1 "BH_BENCHMARK: accretion" "$HERE/run_acc_${name}.log"
done
