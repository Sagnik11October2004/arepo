#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
SOURCE="${SOURCE:-$HERE/output_dm/snap_003.hdf5}"
SYNTH="$HERE/merge_source.hdf5"
PARAM="$HERE/param_merger.txt"
OUT="$HERE/output_merger"
OUTLIST="$HERE/output_list_merger.txt"
NTASKS="${NTASKS:-16}"

if [[ ! -x "$ROOT/ArepoSeedTest" ]]; then
  echo "Missing $ROOT/ArepoSeedTest; run ./build.sh first." >&2
  exit 1
fi
if [[ ! -f "$SOURCE" ]]; then
  echo "Missing one-BH source snapshot: $SOURCE" >&2
  echo "Run the Iteration-10 DM regression first, or set SOURCE=/path/to/snapshot." >&2
  exit 1
fi

python3 "$HERE/make_merger_snapshot.py" "$SOURCE" "$SYNTH"

python3 - "$SYNTH" "$OUTLIST" <<'PY'
import h5py, sys
p, out = sys.argv[1:3]
with h5py.File(p, "r") as f:
    a=float(f["Header"].attrs["Time"])
with open(out, "w") as g:
    g.write(f"{a + 4.0e-6:.15f} 1\n")
    g.write(f"{a + 8.0e-6:.15f} 1\n")
print(f"Merger test starts at a={a:.12f}")
PY

cp "$HERE/param.txt" "$PARAM"
sed -i 's|OutputDir                               ./output|OutputDir                               ./output_merger|' "$PARAM"
sed -i 's|OutputListFilename                      ./output_list.txt|OutputListFilename                      ./output_list_merger.txt|' "$PARAM"
sed -i 's|BHInternalTimestepFactor                1.0|BHInternalTimestepFactor                1.0e30|' "$PARAM"

read TIME_MAX MAX_STEP <<<"$(python3 - "$SYNTH" <<'PY'
import h5py, math, sys
with h5py.File(sys.argv[1],"r") as f:
    a=float(f['Header'].attrs['Time'])
tmax=a+8.0e-6
maxstep=0.1*math.log(tmax/a)
print(f"{tmax:.15f} {maxstep:.15e}")
PY
)"
sed -i "s|TimeMax                                 0.0206000000000000|TimeMax                                 $TIME_MAX|" "$PARAM"
sed -i "s|MaxSizeTimestep                         0.0005|MaxSizeTimestep                         $MAX_STEP|" "$PARAM"
echo "Iteration-11 short restart: TimeMax=$TIME_MAX MaxSizeTimestep=$MAX_STEP"

rm -rf "$OUT"
mkdir -p "$OUT"
cp "$SYNTH" "$OUT/snap_003.hdf5"
# make_merger_snapshot.py already writes the sidecar at
# $HERE/merge_source_expected.json because SYNTH=$HERE/merge_source.hdf5.
# Do not copy the file onto itself under 'set -e'.
EXPECTED="${SYNTH%.hdf5}_expected.json"
if [[ ! -f "$EXPECTED" ]]; then
  echo "Missing merger expectation sidecar: $EXPECTED" >&2
  exit 1
fi
rm -f "$HERE/run_merger.log"

cd "$HERE"
mpirun -np "$NTASKS" "$ROOT/ArepoSeedTest" "$PARAM" 2 3 2>&1 | tee "$HERE/run_merger.log"
