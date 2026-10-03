#!/usr/bin/env python3
"""Verify the Iteration-11 simple BH merger."""

import glob
import json
import math
import os
import re
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "output_merger")
if not os.path.isabs(OUTPUT):
    OUTPUT = os.path.join(HERE, OUTPUT)
LOG = os.path.join(HERE, "run_merger.log")
EXPECTED = os.path.join(HERE, "merge_source_expected.json")

if not os.path.isfile(LOG) or not os.path.isfile(EXPECTED):
    raise SystemExit("FAIL: run merger helper first")

with open(EXPECTED, "r", encoding="utf-8") as fh:
    exp = json.load(fh)
with open(LOG, "r", encoding="utf-8", errors="replace") as fh:
    log = fh.read()

surv_re = re.compile(
    r"BH_FFR: merger survivor ID=(\d+) task=(\d+) members=(\d+) "
    r"Mdyn=([0-9eE+\-.]+) MBH=([0-9eE+\-.]+) Mres=([0-9eE+\-.]+) Mwindbuf=([0-9eE+\-.]+) "
    r"Vx=([0-9eE+\-.]+) Vy=([0-9eE+\-.]+) Vz=([0-9eE+\-.]+)"
)
rows = surv_re.findall(log)
if len(rows) != 1:
    raise SystemExit(f"FAIL: expected exactly one merger survivor event, found {len(rows)}")

pid, task, members, mdyn, mbh, mres, mwbuf, vx, vy, vz = rows[0]
if int(pid) != int(exp["survivor_id"]) or int(members) != 2:
    raise SystemExit("FAIL: deterministic merger survivor/member count mismatch")

if not math.isclose(float(mdyn), exp["merged_mass"], rel_tol=5e-6, abs_tol=1e-12):
    raise SystemExit(f"FAIL: merged dynamical mass {mdyn} != {exp['merged_mass']}")

vlog = np.array([float(vx), float(vy), float(vz)])
vexp = np.asarray(exp["merged_vel"], dtype=float)
if not np.allclose(vlog, vexp, rtol=5e-6, atol=2e-6):
    raise SystemExit(f"FAIL: merger velocity {vlog} != momentum-conserving {vexp}")

consumed_re = re.compile(
    r"BH_FFR: merger consumed ID=(\d+) into ID=(\d+) task=(\d+) "
    r"separationProper=([0-9eE+\-.]+) thresholdProper=([0-9eE+\-.]+)"
)
consumed = consumed_re.findall(log)
if len(consumed) != 1:
    raise SystemExit(f"FAIL: expected exactly one consumed BH, found {len(consumed)}")

loser, survivor, task, sep, threshold = consumed[0]
if int(survivor) != int(exp["survivor_id"]):
    raise SystemExit("FAIL: consumed BH points to wrong survivor")
if not float(sep) < float(threshold):
    raise SystemExit("FAIL: merger fired outside overlap threshold")
if not math.isclose(float(threshold), exp["merge_threshold"], rel_tol=1e-8, abs_tol=1e-12):
    raise SystemExit("FAIL: merger threshold is not 2*BHAccretionRadius")

snaps = sorted(glob.glob(os.path.join(OUTPUT, "snap_*.hdf5")))
post = []
for path in snaps:
    with h5py.File(path, "r") as f:
        a=float(f["Header"].attrs["Time"])
        g=f.get("PartType5")
        n=0 if g is None else len(g["ParticleIDs"])
        print(f"{os.path.basename(path)}: a={a:.9f} N_BH={n}")
        if a > exp["a"] + 1e-10:
            post.append((path,n))

if not post:
    raise SystemExit("FAIL: no post-merger snapshot")
for path,n in post:
    if n != 1:
        raise SystemExit(f"FAIL: {path} contains {n} BHs after merger")
    with h5py.File(path,"r") as f:
        g=f["PartType5"]
        pid=int(g["ParticleIDs"][0])
        mass=float(g["Masses"][0])
        bh=float(g["BH_Mass"][0])
        disk=float(g["BH_DiskMass"][0])
        wbuf=float(g["BH_WindBufferMass"][0])
        if pid != int(exp["survivor_id"]):
            raise SystemExit(f"FAIL: wrong survivor ID {pid}")
        if not math.isclose(mass, bh+disk+wbuf, rel_tol=8e-6, abs_tol=1e-12):
            raise SystemExit("FAIL: post-merger FFR mass ledger does not close")
        if not math.isclose(mass, exp["merged_mass"], rel_tol=8e-6, abs_tol=1e-12):
            raise SystemExit(f"FAIL: post-merger mass {mass} != expected {exp['merged_mass']}")

print()
print("PASS: FFR-MACER Iteration-11 merger test")
print("  two overlapping BHs collapsed to one deterministic survivor")
print("  merger radius is exactly 2 * BHAccretionRadius")
print("  dynamical mass is conserved")
print("  linear momentum is conserved at the merger transaction")
print("  post-merger BHP dynamical-mass ledger closes")
