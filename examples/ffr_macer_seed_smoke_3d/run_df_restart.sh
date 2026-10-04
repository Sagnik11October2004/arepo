#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
SOURCE="${SOURCE:-$HERE/output_dm/snap_003.hdf5}"
SYNTH="$HERE/df_source.hdf5"
PARAM="$HERE/param_df.txt"
OUT="$HERE/output_df"
OUTLIST="$HERE/output_list_df.txt"
NTASKS="${NTASKS:-16}"
EXE="$ROOT/ArepoFFRSnapshotTest"

if [[ ! -x "$EXE" ]]; then
  echo "Missing $EXE; building the snapshot-IC regression executable now." >&2
  "$HERE/build_snapshot_ic.sh"
fi
if [[ ! -x "$EXE" ]]; then
  echo "Build completed without an executable at $EXE." >&2
  exit 1
fi
if [[ ! -f "$SOURCE" ]]; then
  echo "Missing one-BH source snapshot: $SOURCE" >&2
  exit 1
fi

cp "$SOURCE" "$SYNTH"

python3 - "$SYNTH" "$OUTLIST" <<'PY'
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
sed -i 's|InitCondFile                            ./ics|InitCondFile                            ./df_source|' "$PARAM"
sed -i 's|OutputDir                               ./output|OutputDir                               ./output_df|' "$PARAM"
sed -i 's|OutputListFilename                      ./output_list.txt|OutputListFilename                      ./output_list_df.txt|' "$PARAM"
sed -i 's|BHInternalTimestepFactor                1.0|BHInternalTimestepFactor                1.0e30|' "$PARAM"

read TIME_BEGIN TIME_MAX MAX_STEP <<<"$(python3 - "$SYNTH" <<'PY'
import h5py, math, sys
with h5py.File(sys.argv[1],"r") as f:
    a=float(f['Header'].attrs['Time'])
tmax=a+8.0e-6
# Satisfy both the cosmological dln(a) timeline and AREPO's startup check
# MaxSizeTimestep < TimeMax-TimeBegin, which compares the raw parameter values.
dloga=math.log(tmax/a)
da=tmax-a
maxstep=0.1*min(dloga, da)
if not (0.0 < maxstep < da):
    raise SystemExit(f"invalid short-run MaxSizeTimestep={maxstep} for da={da}")
print(f"{a:.15f} {tmax:.15f} {maxstep:.15e}")
PY
)"
sed -i "s|TimeBegin                               0.0200000000000000|TimeBegin                               $TIME_BEGIN|" "$PARAM"
sed -i "s|TimeMax                                 0.0206000000000000|TimeMax                                 $TIME_MAX|" "$PARAM"
sed -i "s|MaxSizeTimestep                         0.0005|MaxSizeTimestep                         $MAX_STEP|" "$PARAM"
echo "Iteration-11 short IC run: TimeBegin=$TIME_BEGIN TimeMax=$TIME_MAX MaxSizeTimestep=$MAX_STEP"

rm -rf "$OUT"
mkdir -p "$OUT"
rm -f "$HERE/run_df.log"

cd "$HERE"
# Start the one-BH hydrodynamic snapshot as a deliberately fresh IC at its
# own scale factor. This uses a dedicated executable without GENERATE_GAS_IN_ICS
# because the source already contains gas. It is not an evolved RestartFlag=2
# continuation and therefore does not bypass the production restart-safety guard.
mpirun -np "$NTASKS" "$EXE" "$PARAM" 2>&1 | tee "$HERE/run_df.log"
