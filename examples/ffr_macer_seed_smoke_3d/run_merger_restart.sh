#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
SOURCE="${SOURCE:-$HERE/output_dm/snap_003.hdf5}"
SYNTH="$HERE/merge_source.hdf5"
SHORT_DA="${SHORT_DA:-2.0e-7}"
SHORT_STEP_FRAC="${SHORT_STEP_FRAC:-0.45}"
PARAM="$HERE/param_merger.txt"
OUT="$HERE/output_merger"
OUTLIST="$HERE/output_list_merger.txt"
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
  echo "Run the Iteration-10 DM regression first, or set SOURCE=/path/to/snapshot." >&2
  exit 1
fi

python3 "$HERE/make_merger_snapshot.py" "$SOURCE" "$SYNTH"

python3 - "$SYNTH" "$OUTLIST" "$SHORT_DA" <<'PY'
import h5py, sys
p, out, da_raw = sys.argv[1:4]
da=float(da_raw)
if not da > 0:
    raise SystemExit(f"SHORT_DA must be positive, got {da}")
with h5py.File(p, "r") as f:
    a=float(f["Header"].attrs["Time"])
with open(out, "w") as g:
    g.write(f"{a + 0.5*da:.15f} 1\n")
    g.write(f"{a + da:.15f} 1\n")
print(f"Merger test starts at a={a:.12f}; SHORT_DA={da:.3e}")
PY

cp "$HERE/param.txt" "$PARAM"
sed -i 's|InitCondFile                            ./ics|InitCondFile                            ./merge_source|' "$PARAM"
sed -i 's|OutputDir                               ./output|OutputDir                               ./output_merger|' "$PARAM"
sed -i 's|OutputListFilename                      ./output_list.txt|OutputListFilename                      ./output_list_merger.txt|' "$PARAM"
sed -i 's|BHInternalTimestepFactor                1.0|BHInternalTimestepFactor                1.0e30|' "$PARAM"
# The synthetic input already contains the required BH(s). Disable high-z FoF
# seeding so this regression does not rebuild a halo catalogue every full step.
sed -i 's|BHSeedMinRedshift                       20.0|BHSeedMinRedshift                       100.0|' "$PARAM"

read TIME_BEGIN TIME_MAX MAX_STEP <<<"$(python3 - "$SYNTH" "$SHORT_DA" "$SHORT_STEP_FRAC" <<'PY'
import h5py, math, sys
with h5py.File(sys.argv[1],"r") as f:
    a=float(f['Header'].attrs['Time'])
da=float(sys.argv[2])
step_frac=float(sys.argv[3])
if not da > 0:
    raise SystemExit(f"SHORT_DA must be positive, got {da}")
if not (0.0 < step_frac < 1.0):
    raise SystemExit(f"SHORT_STEP_FRAC must lie in (0,1), got {step_frac}")
tmax=a+da
# Keep the startup inequality. The merger regression only needs the merge,
# one domain-compaction opportunity, and a post-merger snapshot; it must not
# evolve this synthetic snapshot long enough for unrelated cosmological modes
# to become numerically pathological.
dloga=math.log(tmax/a)
maxstep=step_frac*dloga
if not (0.0 < maxstep < dloga):
    raise SystemExit(f"invalid short-run MaxSizeTimestep={maxstep} for dloga={dloga}")
print(f"{a:.15f} {tmax:.15f} {maxstep:.15e}")
PY
)"
sed -i "s|TimeBegin                               0.0200000000000000|TimeBegin                               $TIME_BEGIN|" "$PARAM"
sed -i "s|TimeMax                                 0.0206000000000000|TimeMax                                 $TIME_MAX|" "$PARAM"
sed -i "s|MaxSizeTimestep                         0.0005|MaxSizeTimestep                         $MAX_STEP|" "$PARAM"
echo "Iteration-11 fast short IC run: TimeBegin=$TIME_BEGIN TimeMax=$TIME_MAX MaxSizeTimestep=$MAX_STEP SHORT_DA=$SHORT_DA STEP_FRAC=$SHORT_STEP_FRAC"

rm -rf "$OUT"
mkdir -p "$OUT"
# make_merger_snapshot.py already writes the sidecar at
# $HERE/merge_source_expected.json because SYNTH=$HERE/merge_source.hdf5.
EXPECTED="${SYNTH%.hdf5}_expected.json"
if [[ ! -f "$EXPECTED" ]]; then
  echo "Missing merger expectation sidecar: $EXPECTED" >&2
  exit 1
fi
rm -f "$HERE/run_merger.log"

cd "$HERE"
# Deliberately start the synthetic two-BH hydrodynamic file as fresh IC state
# with the dedicated snapshot-IC executable (no GENERATE_GAS_IN_ICS). The
# production RestartFlag=2 guard remains intact for true evolved continuations.
mpirun -np "$NTASKS" "$EXE" "$PARAM" 2>&1 | tee "$HERE/run_merger.log"