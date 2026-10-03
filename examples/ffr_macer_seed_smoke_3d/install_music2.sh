#!/usr/bin/env bash
set -euo pipefail

DEST="${MUSIC2_DIR:-$HOME/Research/MUSIC2}"
JOBS="${JOBS:-$(nproc)}"

if [[ ! -d "$DEST/.git" ]]; then
  git clone https://github.com/cosmo-sims/MUSIC2.git "$DEST"
else
  git -C "$DEST" pull --ff-only
fi

cmake -S "$DEST" -B "$DEST/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$DEST/build" -j"$JOBS"

echo
echo "MUSIC2 executable expected at: $DEST/build/MUSIC"
