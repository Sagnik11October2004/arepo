#!/usr/bin/env python3
"""Verify the Iteration-11 Chandrasekhar dynamical-friction kick."""

import math
import os
import re

HERE=os.path.dirname(os.path.abspath(__file__))
LOG=os.path.join(HERE,"run_df.log")
if not os.path.isfile(LOG):
    raise SystemExit("FAIL: run_df.log not found; run ./run_df_restart.sh first")

with open(LOG,"r",encoding="utf-8",errors="replace") as fh:
    log=fh.read()

pat=re.compile(
    r"BH_FFR: dynamical friction ID=(\d+) task=(\d+) "
    r"rhoDM=([0-9eE+\-.]+) sigmaDM=([0-9eE+\-.]+) "
    r"vrel0=([0-9eE+\-.]+) vrel1=([0-9eE+\-.]+) "
    r"tdfCode=([0-9eE+\-.]+) dtCode=([0-9eE+\-.]+) "
    r"frac=([0-9eE+\-.]+) lnLambda=([0-9eE+\-.]+)"
)
rows=pat.findall(log)
if not rows:
    raise SystemExit("FAIL: no dynamical-friction kick was logged")

strict=0
for row in rows:
    pid=row[0]
    rho,sigma,v0,v1,tdf,dt,frac,lnl=map(float,row[2:])
    if not all(math.isfinite(x) for x in (rho,sigma,v0,v1,tdf,dt,frac,lnl)):
        raise SystemExit(f"FAIL: BH {pid} non-finite DF diagnostic")
    if not (rho>0 and sigma>0 and v0>0 and tdf>0 and dt>0):
        raise SystemExit(f"FAIL: BH {pid} non-positive DF environment/timescale")
    if not (0 <= frac <= 0.5 + 1e-12):
        raise SystemExit(f"FAIL: BH {pid} damping fraction {frac} outside [0,0.5]")
    if v1 > v0*(1+2e-10):
        raise SystemExit(f"FAIL: BH {pid} dynamical friction increased vrel {v0}->{v1}")
    if frac>0 and v1 < v0:
        strict += 1
    expected=v0*(1-frac)
    if not math.isclose(v1,expected,rel_tol=3e-5,abs_tol=2e-7):
        raise SystemExit(f"FAIL: BH {pid} vrel1={v1} != capped drag expectation {expected}")
    if not math.isclose(lnl,3.0,rel_tol=0,abs_tol=1e-12):
        raise SystemExit(f"FAIL: BH {pid} Coulomb log {lnl} != 3")

if strict==0:
    raise SystemExit("FAIL: DF module never strictly reduced BH velocity relative to local DM")

if "BH_FFR: merger survivor" in log:
    raise SystemExit("FAIL: merger occurred in the single-BH DF regression")

print("PASS: FFR-MACER Iteration-11 dynamical-friction test")
print(f"  DF kicks checked             : {len(rows)}")
print(f"  strict relative-speed drops : {strict}")
print("  Chandrasekhar drag always opposes BH motion relative to local DM")
print("  per-update damping is capped at 50 percent")
print("  no merger contaminated the single-BH regression")
