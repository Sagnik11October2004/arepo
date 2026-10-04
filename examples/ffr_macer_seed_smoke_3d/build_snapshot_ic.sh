#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
CONFIG="examples/ffr_macer_seed_smoke_3d/Config_snapshot_ic.sh"
BUILD_DIR="build_ffr_snapshot"
EXEC="ArepoFFRSnapshotTest"
JOBS="${JOBS:-$(nproc)}"

cd "$ROOT"
make clean CONFIG="$CONFIG" BUILD_DIR="$BUILD_DIR" EXEC="$EXEC"
make -j"$JOBS" CONFIG="$CONFIG" BUILD_DIR="$BUILD_DIR" EXEC="$EXEC"

echo
echo "Built: $ROOT/$EXEC"
