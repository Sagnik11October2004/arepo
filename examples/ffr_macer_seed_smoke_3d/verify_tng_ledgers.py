#!/usr/bin/env python3
import argparse
import math
import re
from pathlib import Path

KV=re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=([^\s]+)")

def main():
    ap=argparse.ArgumentParser(description="Verify cumulative TNG model-energy accounting from a run log.")
    ap.add_argument("log")
    ap.add_argument("--tol", type=float, default=5e-10)
    args=ap.parse_args()
    text=Path(args.log).read_text(errors="replace")
    rows=[]
    for line in text.splitlines():
        if "BH_TNG: coupling-ledger" not in line:
            continue
        d=dict(KV.findall(line))
        rel=float(d["rel"])
        active=float(d["activeFrac"])
        if not math.isfinite(rel) or abs(rel) > args.tol:
            raise AssertionError(f"energy closure failed: {line}")
        if not (0.0 <= active <= 1.0 + 1e-12):
            raise AssertionError(f"invalid active fraction: {line}")
        rows.append((active,int(d["thermAgeTransactions"]),int(d["kinAgeTransactions"]),rel))
    if not rows:
        raise SystemExit("No BH_TNG: coupling-ledger diagnostics found")
    print(f"N={len(rows)} max|energy closure rel|={max(abs(x[3]) for x in rows):.3e} "
          f"min activeFrac={min(x[0] for x in rows):.6g} "
          f"max thermal wait={max(x[1] for x in rows)} transactions "
          f"max kinetic wait={max(x[2] for x in rows)} transactions")
    print("TNG cumulative energy ledger verification: PASS")

if __name__ == "__main__":
    main()
