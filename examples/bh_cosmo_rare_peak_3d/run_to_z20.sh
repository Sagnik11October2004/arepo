#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
NTASKS="${NTASKS:-8}"
MPIEXEC="${MPIEXEC:-mpirun}"

cd "$HERE"

MIN_FREE_GIB="${MIN_FREE_GIB:-5}"
FREE_KIB="$(df -Pk "$HERE" | awk 'NR==2 {print $4}')"
FREE_GIB="$(awk -v k="$FREE_KIB" 'BEGIN {printf "%.2f", k/1024/1024}')"
echo "Filesystem free space at $HERE: ${FREE_GIB} GiB"
awk -v k="$FREE_KIB" -v min="$MIN_FREE_GIB" 'BEGIN {exit !(k >= min*1024*1024)}' || {
  echo "ERROR: less than ${MIN_FREE_GIB} GiB free; refusing to start the scout." >&2
  exit 1
}

if [[ -s ics.hdf5 ]]; then
  if ! python3 verify_ic.py ics.hdf5; then
    echo "Existing ics.hdf5 does not match the z=49 64^3 scout; regenerating it."
    ./generate_ics.sh
  fi
else
  echo "No ics.hdf5 found; generating z=49 64^3 MUSIC scout IC."
  ./generate_ics.sh
fi

if [[ ! -x "$ROOT/ArepoRarePeakInitial" || ! -x "$ROOT/ArepoRarePeak" ]]; then
  echo "Rare-peak executables missing; building them."
  ./build.sh
fi

mkdir -p output_preseed output_seed_window logs ics

if [[ "${RESUME_PRESEED:-0}" == "1" ]]; then
  echo
  echo "=== Stage A resume: latest existing BH-free snapshot -> z=22 ==="

  LATEST_PRESEED="$(ls output_preseed/snap_*.hdf5 2>/dev/null | sort | tail -n 1 || true)"
  if [[ -z "$LATEST_PRESEED" ]]; then
    echo "ERROR: RESUME_PRESEED=1 but no output_preseed/snap_*.hdf5 exists." >&2
    exit 1
  fi

  SNAPNUM="$(basename "$LATEST_PRESEED" .hdf5 | sed 's/^snap_//')"
  SNAPNUM_DEC="$((10#$SNAPNUM))"

  python3 - "$LATEST_PRESEED" <<'PY'
import h5py, sys
p=sys.argv[1]
with h5py.File(p,"r") as f:
    z=float(f["Header"].attrs["Redshift"])
    nbh=len(f["PartType5"]["ParticleIDs"]) if "PartType5" in f else 0
if not (22.0 - 2e-3 <= z < 99.0):
    raise SystemExit(f"Refusing Stage-A resume from unexpected redshift: {p} z={z}")
if nbh != 0:
    raise SystemExit(f"Stage-A resume snapshot must be BH-free; found {nbh} BH(s)")
print(f"Resuming Stage A from {p}: z={z:.6f}, N_BH={nbh}")
PY

  rm -rf output_seed_window
  mkdir -p output_seed_window
  : > logs/preseed_resume.log

  "$MPIEXEC" -np "$NTASKS" "$ROOT/ArepoRarePeak" param_preseed.txt 2 "$SNAPNUM_DEC" \
    </dev/null 2>&1 | tee logs/preseed_resume.log
else
  EXISTING_PRESEED="$(ls output_preseed/snap_*.hdf5 2>/dev/null | sort | tail -n 1 || true)"
  if [[ -n "$EXISTING_PRESEED" && "${FORCE_FRESH_PRESEED:-0}" != "1" ]]; then
    echo "ERROR: existing Stage-A snapshot found: $EXISTING_PRESEED" >&2
    echo "Refusing to delete it." >&2
    echo "Use RESUME_PRESEED=1 to continue, or FORCE_FRESH_PRESEED=1 only if you deliberately want a fresh z=49 scout run." >&2
    exit 1
  fi

  rm -rf output_preseed output_seed_window logs
  mkdir -p output_preseed output_seed_window logs ics

  echo
  echo "=== Stage A: quick-scout BH-free cosmology, z=49 -> z=22 ==="
  "$MPIEXEC" -np "$NTASKS" "$ROOT/ArepoRarePeakInitial" param_preseed.txt \
    </dev/null 2>&1 | tee logs/preseed.log
fi

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

python3 audit_seed_event.py \
  logs/seed_window.log \
  ics/common_seed_metadata.json \
  --z-min 20 --z-max 22 --expected-seed-mass 1.0e5

Z20_SNAP="$(ls output_seed_window/snap_*.hdf5 2>/dev/null | sort | tail -n 1)"
if [[ -z "$Z20_SNAP" ]]; then
  echo "ERROR: no z=20 snapshot found." >&2
  exit 1
fi

python3 - "$Z20_SNAP" ics/common_seed_metadata.json <<'PY'
import json, h5py, numpy as np, sys
p=sys.argv[1]
meta_path=sys.argv[2]
h=0.68
meta=json.load(open(meta_path))
with h5py.File(p,"r") as f:
    z=float(f["Header"].attrs["Redshift"])
    if "PartType5" not in f:
        raise SystemExit("No BH was seeded in 20<z<22. Inspect logs/seed_window.log and the 1e7 Msun host threshold.")
    g=f["PartType5"]
    n=len(g["ParticleIDs"])
    ids=np.asarray(g["ParticleIDs"], dtype=np.uint64)
    masses=np.asarray(g["Masses"], dtype=float)*1e10/h
if abs(z-20.0)>2e-3:
    raise SystemExit(f"Latest seed-window snapshot is not z=20: z={z}")
if n != 1:
    raise SystemExit(f"Single-seed run expected exactly one BH; found {n}")
if int(ids[0]) != int(meta["bh_id"]):
    raise SystemExit(f"Final BH ID {int(ids[0])} does not match seed-event ID {meta['bh_id']}")
if abs(masses[0]/1e5-1.0)>2e-6:
    raise SystemExit(f"Capture-off seed should remain 1e5 Msun; found {masses[0]} Msun")
print(
    f"Validated common seeded state: {p} z={z:.6f} N_BH={n} "
    f"M_BH={masses[0]:.6g} Msun BH_ID={int(ids[0])} "
    f"seed_z={meta['seed_redshift']:.6f} "
    f"Mhost={meta['host_fof_mass_msun']:.6e} Msun"
)
PY

python3 make_common_branch_ic.py "$Z20_SNAP" ics/common_z20_seeded.hdf5

echo
echo "Rare-peak common state is ready:"
echo "  $HERE/ics/common_z20_seeded.hdf5"
echo "  $HERE/ics/common_seed_metadata.json"
echo "Use prepare_branch_params.py / run_branches.sh for post-z=20 science branches."
