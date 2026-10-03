#!/usr/bin/env python3
"""Verify Iteration-10 cached full-gravity-tree SigmaDM."""

import glob
import math
import os
import re
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "output_dm")
if not os.path.isabs(OUTPUT):
    OUTPUT = os.path.join(HERE, OUTPUT)
LOG = os.path.join(HERE, "run_dm.log")

REQUESTED = 64
FW = 0.01
FJ = 0.10


def fail(msg):
    raise SystemExit("FAIL: " + msg)


if not os.path.isfile(LOG):
    fail("run_dm.log not found; run ./run_dm_restart.sh first")

with open(LOG, "r", encoding="utf-8", errors="replace") as fh:
    log = fh.read()

cache_re = re.compile(
    r"BH_FFR: sigmaDM cache ID=(\d+) task=(\d+) source=([^\s]+) "
    r"count=(\d+) requested=(\d+) sigma=([0-9eE+\-.]+) rmaxProper=([0-9eE+\-.]+)"
)
cache_rows = cache_re.findall(log)
if not cache_rows:
    fail("no SigmaDM cache update was logged")

dm_rows = []
for pid, task, source, ctxt, rtxt, stxt, rmaxtxt in cache_rows:
    count = int(ctxt)
    requested = int(rtxt)
    sigma = float(stxt)
    rmax = float(rmaxtxt)

    if requested != REQUESTED:
        fail(f"BH {pid}: requested DM neighbours {requested}, expected {REQUESTED}")
    if source == "DM":
        dm_rows.append((pid, sigma))
        if count != REQUESTED:
            fail(f"BH {pid}: exact DM sample has count={count}, expected {REQUESTED}")
        if not (math.isfinite(sigma) and sigma > 0):
            fail(f"BH {pid}: invalid positive DM dispersion sigma={sigma}")
        if not (math.isfinite(rmax) and rmax > 0):
            fail(f"BH {pid}: invalid nearest-DM sample radius {rmax}")
    elif source == "gas-fallback":
        fail(f"BH {pid}: cosmological smoke test unexpectedly used gas fallback")
    else:
        fail(f"BH {pid}: unknown SigmaDM source {source}")

if not dm_rows:
    fail("no successful DM-sourced SigmaDM cache update")

threshold_re = re.compile(
    r"BH_FFR: binding threshold ID=(\d+) task=(\d+) "
    r"Menc=([0-9eE+\-.]+) sigmaDM=([0-9eE+\-.]+) "
    r"vbind2=([0-9eE+\-.]+) central=(\d+) "
    r"EthWind=([0-9eE+\-.]+) EthJet=([0-9eE+\-.]+)"
)
threshold_rows = threshold_re.findall(log)
if not threshold_rows:
    fail("no positive SigmaDM binding threshold was logged")

checked_thresholds = 0
for pid, task, mtxt, stxt, vtxt, central_txt, ewtxt, ejtxt in threshold_rows:
    menc = float(mtxt)
    sigma = float(stxt)
    vbind2 = float(vtxt)
    central = int(central_txt)
    ethw = float(ewtxt)
    ethj = float(ejtxt)

    if central != 0:
        fail(f"BH {pid}: DM-only regression has central binding term enabled")
    if not all(math.isfinite(x) and x > 0 for x in (menc, sigma, vbind2, ethw, ethj)):
        fail(f"BH {pid}: invalid positive threshold diagnostics")

    if not math.isclose(vbind2, sigma * sigma, rel_tol=4e-5, abs_tol=1e-12):
        fail(f"BH {pid}: vbind2={vbind2} does not equal sigmaDM^2={sigma*sigma}")

    base = 0.5 * menc * vbind2
    if not math.isclose(ethw, FW * base, rel_tol=6e-5, abs_tol=1e-14):
        fail(f"BH {pid}: EthWind={ethw} inconsistent with DM-only binding scale {FW*base}")
    if not math.isclose(ethj, FJ * base, rel_tol=6e-5, abs_tol=1e-14):
        fail(f"BH {pid}: EthJet={ethj} inconsistent with DM-only binding scale {FJ*base}")
    checked_thresholds += 1

if "BH_FFR: wind packet fired" in log or "BH_FFR: jet packet fired" in log:
    fail("mechanical feedback fired even though capture/inner power are disabled")

snaps = sorted(glob.glob(os.path.join(OUTPUT, "snap_*.hdf5")))
if not snaps:
    fail(f"no snapshots found in {OUTPUT}")

seen_bh = False
seen_positive_sigma = False
last_sigma_by_id = {}

for path in snaps:
    with h5py.File(path, "r") as f:
        a = float(f["Header"].attrs["Time"])
        z = float(f["Header"].attrs["Redshift"])
        g = f.get("PartType5")
        nbh = 0 if g is None else len(g["ParticleIDs"])
        print(f"{os.path.basename(path)}: a={a:.9f} z={z:.6f} N_BH={nbh}")
        if nbh == 0:
            continue

        seen_bh = True
        for name in (
            "BH_SigmaDM", "BH_EthWind", "BH_EthJet",
            "BH_MdotSupply", "BH_MdotFeed", "BH_MdotHorizon", "BH_MdotWind",
            "BH_EWind", "BH_EJet",
        ):
            if name not in g:
                fail(f"{os.path.basename(path)} missing {name}")

        ids = np.asarray(g["ParticleIDs"][:], dtype=np.uint64)
        for i, pid0 in enumerate(ids):
            pid = int(pid0)
            sigma = float(g["BH_SigmaDM"][i])
            ethw = float(g["BH_EthWind"][i])
            ethj = float(g["BH_EthJet"][i])

            if not all(math.isfinite(x) and x >= 0 for x in (sigma, ethw, ethj)):
                fail(f"{os.path.basename(path)} BH {pid}: invalid DM diagnostics")

            if sigma > 0:
                seen_positive_sigma = True
                if not (ethw > 0 and ethj > 0):
                    fail(f"{os.path.basename(path)} BH {pid}: positive SigmaDM but non-positive thresholds")
                if not math.isclose(ethj / ethw, FJ / FW, rel_tol=2e-5):
                    fail(f"{os.path.basename(path)} BH {pid}: threshold factor ratio mismatch")

            for name in ("BH_MdotSupply", "BH_MdotFeed", "BH_MdotHorizon", "BH_MdotWind", "BH_EWind", "BH_EJet"):
                val = float(g[name][i])
                if abs(val) > 1e-12:
                    fail(f"{os.path.basename(path)} BH {pid}: {name}={val} should remain zero in DM-only regression")

            last_sigma_by_id[pid] = sigma
            print(
                f"  BH {pid}: SigmaDM={sigma:.8e} "
                f"EthWind={ethw:.8e} EthJet={ethj:.8e}"
            )

if not seen_bh:
    fail("no BH formed")
if not seen_positive_sigma:
    fail("BH formed but cached SigmaDM never became positive")

# The snapshot cache should agree with one of the logged full-tree values for
# that BH to log precision. It need not be the final log line if a later
# full-tree update occurs after the final snapshot write.
for pid, sigma_snap in last_sigma_by_id.items():
    logged = [sigma for row_pid, sigma in dm_rows if int(row_pid) == pid]
    if sigma_snap > 0 and not any(math.isclose(sigma_snap, x, rel_tol=2e-5, abs_tol=1e-8) for x in logged):
        fail(f"BH {pid}: snapshot SigmaDM={sigma_snap} does not match any logged cache value")

print()
print("PASS: FFR-MACER Iteration-10 DM-dispersion test")
print(f"  full-tree DM cache updates       : {len(dm_rows)}")
print(f"  DM-only threshold checks         : {checked_thresholds}")
print(f"  exact nearest-DM sample size     : {REQUESTED}")
print("  cached SigmaDM is finite and positive")
print("  no gas fallback was needed in the cosmological smoke box")
print("  with central binding disabled, vbind^2 = SigmaDM^2")
print("  wind/jet burst thresholds match the DM binding-energy scale")
print("  no capture or mechanical feedback contaminated the regression")
