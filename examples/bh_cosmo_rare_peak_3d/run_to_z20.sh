#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
NTASKS="${NTASKS:-8}"
MPIEXEC="${MPIEXEC:-mpirun}"

cd "$HERE"

if [[ ! -s ics.hdf5 ]]; then
  echo "No ics.hdf5 found; generating z=99 MUSIC IC."
  ./generate_ics.sh
fi

if [[ ! -x "$ROOT/ArepoRarePeakInitial" || ! -x "$ROOT/ArepoRarePeak" ]]; then
  echo "Rare-peak executables missing; building them."
  ./build.sh
fi

rm -rf output_preseed output_seed_window logs
mkdir -p output_preseed output_seed_window logs ics

echo
echo "=== Stage A: common BH-free cosmology, z=99 -> z=22 ==="
"$MPIEXEC" -np "$NTASKS" "$ROOT/ArepoRarePeakInitial" param_preseed.txt \
  </dev/null 2>&1 | tee logs/preseed.log

if [[ ! -f output_preseed/end ]]; then
  echo "ERROR: pre-seed stage did not finish cleanly." >&2
  exit 1
fi

Z22_SNAP="$(ls output_preseed/snap_*.hdf5 2>/dev/null | sort | tail -n 1)"
if [[ -z "$Z22_SNAP" ]]; then
  echo "ERROR: no z=22 snapshot found." >&2
  exit 1
fi

python3 - "$Z22_SNAP" <<'PY'
import h5py, sys
p=sys.argv[1]
with h5py.File(p,"r") as f:
    z=float(f["Header"].attrs["Redshift"])
    nbh=len(f["PartType5"]["ParticleIDs"]) if "PartType5" in f else 0
if abs(z-22.0)>2e-3:
    raise SystemExit(f"Latest pre-seed snapshot is not z=22: {p} z={z}")
if nbh != 0:
    raise SystemExit(f"Pre-seed checkpoint must be BH-free; found {nbh} BH(s)")
print(f"Validated BH-free z=22 checkpoint: {p} z={z:.6f}")
PY

cp -f "$Z22_SNAP" ics/z22_preseed.hdf5

echo
echo "=== Stage B: single-most-massive FoF seeding window, z=22 -> z=20 ==="
echo "    Capture and feedback remain OFF; only the seed transaction is active."
"$MPIEXEC" -np "$NTASKS" "$ROOT/ArepoRarePeak" param_seed_window.txt \
  </dev/null 2>&1 | tee logs/seed_window.log

if [[ ! -f output_seed_window/end ]]; then
  echo "ERROR: seed-window stage did not finish cleanly." >&2
  exit 1
fi

Z20_SNAP="$(ls output_seed_window/snap_*.hdf5 2>/dev/null | sort | tail -n 1)"
if [[ -z "$Z20_SNAP" ]]; then
  echo "ERROR: no z=20 snapshot found." >&2
  exit 1
fi

python3 - "$Z20_SNAP" <<'PY'
import h5py, numpy as np, sys
p=sys.argv[1]
h=0.68
with h5py.File(p,"r") as f:
    z=float(f["Header"].attrs["Redshift"])
    if "PartType5" not in f:
        raise SystemExit("No BH was seeded in 20<z<22. Inspect logs/seed_window.log and the 1e7 Msun host threshold.")
    g=f["PartType5"]
    n=len(g["ParticleIDs"])
    masses=np.asarray(g["Masses"], dtype=float)*1e10/h
if abs(z-20.0)>2e-3:
    raise SystemExit(f"Latest seed-window snapshot is not z=20: z={z}")
if n != 1:
    raise SystemExit(f"Single-seed run expected exactly one BH; found {n}")
if abs(masses[0]/1e5-1.0)>2e-6:
    raise SystemExit(f"Capture-off seed should remain 1e5 Msun; found {masses[0]} Msun")
print(f"Validated common seeded state: {p} z={z:.6f} N_BH={n} M_BH={masses[0]:.6g} Msun")
PY

python3 make_common_branch_ic.py "$Z20_SNAP" ics/common_z20_seeded.hdf5

echo
echo "Rare-peak common state is ready:"
echo "  $HERE/ics/common_z20_seeded.hdf5"
echo "Use prepare_branch_params.py / run_branches.sh for post-z=20 science branches."
