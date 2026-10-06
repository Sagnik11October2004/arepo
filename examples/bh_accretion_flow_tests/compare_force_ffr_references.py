#!/usr/bin/env python3
"""Historical entry point; forwards to the current comparison with references."""
import subprocess
import sys
from pathlib import Path
here = Path(__file__).resolve().parent
raise SystemExit(subprocess.call([
    sys.executable, str(here / "compare_flow_models.py"), "--with-reference", *sys.argv[1:]
]))
