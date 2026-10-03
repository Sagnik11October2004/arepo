#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
SOURCE="${SOURCE:-$HERE/output_dm/snap_003.hdf5}"
PARAM="$HERE/param_df.txt"
OUT="$HERE/output_df"
OUTLIST="$HERE/output_list_df.txt"
NTASKS="${NTASKS:-16}"

if [[ ! -x "$ROOT/ArepoSeedTest" ]]; then
  echo "Missing $ROOT/ArepoSeedTest; run ./build.sh first." >&2
  exit 1
fi
if [[ ! -f "$SOURCE" ]]; then
  echo "Missing one-BH source snapshot: $SOURCE" >&2
  exit 1
fi

python3 - "$SOURCE" "$OUTLIST" <<'PY'
import h5py,sys
p,out=sys.argv[1:3]
with h5py.File(p,"r") as f:
    a=float(f["Header"].attrs["Time"])
    n=len(f["PartType5"]["ParticleIDs"]) if "PartType5" in f else 0
if n != 1:
    raise SystemExit(f"need exactly one BH, found {n}")
with open(out,"w") as g:
    g.write(f"{a + 4.0e-6:.15f} 1\n")
    g.write(f"{a + 8.0e-6:.15f} 1\n")
print(f"DF test starts at a={a:.12f}")
PY

cp "$HERE/param.txt" "$PARAM"
sed -i 's|OutputDir                               ./output|OutputDir                               ./output_df|' "$PARAM"
sed -i 's|OutputListFilename                      ./output_list.txt|OutputListFilename                      ./output_list_df.txt|' "$PARAM"
sed -i 's|BHInternalTimestepFactor                1.0|BHInternalTimestepFactor                1.0e30|' "$PARAM"

read TIME_MAX MAX_STEP <<<"$(python3 - "$SOURCE" <<'PY'
import h5py, math, sys
with h5py.File(sys.argv[1],"r") as f:
    a=float(f['Header'].attrs['Time'])
tmax=a+8.0e-6
# RestartFlag=2 resets TimeBegin to a. Keep the test timestep comfortably
# below the resulting full logarithmic timeline span, rather than inheriting
# the canonical 5e-4 dloga step from the much longer smoke run.
maxstep=0.1*math.log(tmax/a)
print(f"{tmax:.15f} {maxstep:.15e}")
PY
)"
sed -i "s|TimeMax                                 0.0206000000000000|TimeMax                                 $TIME_MAX|" "$PARAM"
sed -i "s|MaxSizeTimestep                         0.0005|MaxSizeTimestep                         $MAX_STEP|" "$PARAM"
echo "Iteration-11 short restart: TimeMax=$TIME_MAX MaxSizeTimestep=$MAX_STEP"

rm -rf "$OUT"
mkdir -p "$OUT"
cp "$SOURCE" "$OUT/snap_003.hdf5"
rm -f "$HERE/run_df.log"

cd "$HERE"
mpirun -np "$NTASKS" "$ROOT/ArepoSeedTest" "$PARAM" 2 3 2>&1 | tee "$HERE/run_df.log"
