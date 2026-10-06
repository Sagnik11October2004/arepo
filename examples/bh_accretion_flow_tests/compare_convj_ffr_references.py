#!/usr/bin/env python3
"""Reference-oriented compatibility entry point for current mode 5."""
import subprocess
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
cmd = [sys.executable, str(here / "compare_flow_models.py"), "--with-reference", *sys.argv[1:]]
raise SystemExit(subprocess.call(cmd))
