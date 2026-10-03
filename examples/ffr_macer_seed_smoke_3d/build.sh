#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
CONFIG="examples/ffr_macer_seed_smoke_3d/Config.sh"
JOBS="${JOBS:-$(nproc)}"

cd "$ROOT"
make clean CONFIG="$CONFIG" EXEC=ArepoSeedTest
make -j"$JOBS" CONFIG="$CONFIG" EXEC=ArepoSeedTest

echo
echo "Built: $ROOT/ArepoSeedTest"
