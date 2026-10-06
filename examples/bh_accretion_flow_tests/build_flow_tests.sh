#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
CONFIG="examples/bh_accretion_flow_tests/Config_flow_tests.sh"
BUILD_DIR="build_bh_flows"
EXEC="ArepoBHFlows"
JOBS="${JOBS:-$(nproc)}"

cd "$ROOT"
rm -rf "$BUILD_DIR"
rm -f "$EXEC"

make -j"$JOBS" CONFIG="$CONFIG" BUILD_DIR="$BUILD_DIR" EXEC="$EXEC"

echo
echo "Built: $ROOT/$EXEC"
grep -E '^(DOUBLEPRECISION|DOUBLEPRECISION_FFTW|INPUT_IN_DOUBLEPRECISION|OUTPUT_IN_DOUBLEPRECISION|OUTPUT_COORDINATES_IN_DOUBLEPRECISION|NGB_TREE_DOUBLEPRECISION|HAVE_HDF5|BLACKHOLE_FFR|FOF)' "$ROOT/$CONFIG"
