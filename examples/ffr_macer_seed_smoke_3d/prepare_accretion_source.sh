#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ $# -ne 1 ]]; then
  echo "Usage: $0 /path/to/one_bh_source_snapshot.hdf5" >&2
  exit 2
fi

SRC="$(python3 - "$1" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).resolve())
PY
)"
DST="$HERE/accretion_source.hdf5"
TMP="$HERE/.accretion_source.hdf5.tmp"

if [[ ! -f "$SRC" ]]; then
  echo "Missing source snapshot: $SRC" >&2
  exit 1
fi
case "$SRC" in
  *.hdf5) ;;
  *) echo "Source must be an HDF5 snapshot: $SRC" >&2; exit 1 ;;
esac

python3 - "$SRC" <<'PY'
import h5py, sys
p=sys.argv[1]
with h5py.File(p,"r") as f:
    if "PartType5" not in f:
        raise SystemExit(f"{p}: missing PartType5")
    n=len(f["PartType5"]["ParticleIDs"])
    if n != 1:
        raise SystemExit(f"{p}: expected exactly one BH, found {n}")
    disk=f["PartType5"].get("BH_DiskMass")
    wind=f["PartType5"].get("BH_WindBufferMass")
    if disk is not None and float(disk[0]) != 0.0:
        raise SystemExit(f"{p}: BH_DiskMass must be zero for direct-backend source")
    if wind is not None and float(wind[0]) != 0.0:
        raise SystemExit(f"{p}: BH_WindBufferMass must be zero for direct-backend source")
PY

rm -f "$TMP"
cp --reflink=auto --sparse=always "$SRC" "$TMP"
mv -f "$TMP" "$DST"

if [[ -L "$DST" ]]; then
  echo "Internal error: durable destination became a symlink." >&2
  exit 1
fi

echo "Prepared standalone benchmark source:"
echo "  $DST"
ls -lh "$DST"
