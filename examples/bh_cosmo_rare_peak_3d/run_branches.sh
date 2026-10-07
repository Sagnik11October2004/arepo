#!/usr/bin/env bash
set -euo pipefail

cat >&2 <<'EOF'
This legacy snapshot-based branch runner is disabled.

BH-containing evolved FFR/MACER states must be branched from the fully
serialized native RestartFlag=1 checkpoint, not from a snapshot/IC restart.
The old runner also deleted output directories before launch.

Use instead:

  BRANCHES="convj_shell_ffr__macer" NTASKS=4 ./run_native_restart_suite.sh

See README.md and params_native/ for the checked-in science branch suite.
EOF
exit 2
