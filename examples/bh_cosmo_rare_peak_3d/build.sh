#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
JOBS="${JOBS:-$(nproc)}"

cd "$ROOT"

make clean CONFIG=examples/bh_cosmo_rare_peak_3d/Config_initial.sh \
  BUILD_DIR=build_rarepeak_initial EXEC=ArepoRarePeakInitial
make -j"$JOBS" CONFIG=examples/bh_cosmo_rare_peak_3d/Config_initial.sh \
  BUILD_DIR=build_rarepeak_initial EXEC=ArepoRarePeakInitial

make clean CONFIG=examples/bh_cosmo_rare_peak_3d/Config_evolved.sh \
  BUILD_DIR=build_rarepeak_evolved EXEC=ArepoRarePeak
make -j"$JOBS" CONFIG=examples/bh_cosmo_rare_peak_3d/Config_evolved.sh \
  BUILD_DIR=build_rarepeak_evolved EXEC=ArepoRarePeak

echo
echo "Built:"
echo "  $ROOT/ArepoRarePeakInitial  (DM-only MUSIC IC -> gas+DM split)"
echo "  $ROOT/ArepoRarePeak         (already-split cosmological snapshots/branch ICs)"
