#!/usr/bin/env python3
"""Verify Iteration-7 bursty bipolar wind feedback."""

import glob
import math
import os
import re
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "output_feedback")
if not os.path.isabs(OUTPUT):
    OUTPUT = os.path.join(HERE, OUTPUT)
LOG = os.path.join(HERE, "run_feedback.log")

HUBBLE = 0.68
UNIT_MASS_MSUN = 1.0e10
UNIT_LENGTH_CM = 3.085678e21
UNIT_VELOCITY_CM_S = 1.0e5
SEC_PER_YEAR = 365.25 * 24.0 * 3600.0
UNIT_TIME_S = UNIT_LENGTH_CM / UNIT_VELOCITY_CM_S

MASS_TO_MSUN = UNIT_MASS_MSUN / HUBBLE
RATE_TO_MSUN_YR = UNIT_MASS_MSUN / (UNIT_TIME_S / SEC_PER_YEAR)

WIND_FACTOR = 0.001
JET_FACTOR = 1.0e30
MIN_TARGETS = 8


def fail(msg):
    raise SystemExit("FAIL: " + msg)


if not os.path.isfile(LOG):
    fail("run_feedback.log not found; run ./run_feedback_restart.sh first")

with open(LOG, "r", encoding="utf-8", errors="replace") as fh:
    log = fh.read()

event_re = re.compile(
    r"BH_FFR: wind packet fired ID=(\d+) task=(\d+) "
    r"E=([0-9eE+\-.]+) erg dM=([0-9eE+\-.]+) Msun "
    r"q=([0-9eE+\-.]+) Nplus=(\d+) Nminus=(\d+) "
    r"dErel=([0-9eE+\-.]+) pbal=([0-9eE+\-.]+)"
)
events = event_re.findall(log)
if not events:
    fail("no resolved wind packet was fired")
if "BH_FFR: jet packet fired" in log:
    fail("jet packet fired during the isolated Iteration-7 wind regression")

for pid, task, etxt, mtxt, qtxt, nptxt, nmtxt, detxt, pbtxt in events:
    e = float(etxt)
    dm = float(mtxt)
    q = float(qtxt)
    np_ = int(nptxt)
    nm_ = int(nmtxt)
    de = float(detxt)
    pb = float(pbtxt)

    if not (math.isfinite(e) and e > 0):
        fail(f"BH {pid}: non-positive packet energy {e}")
    if not (math.isfinite(dm) and dm >= 0):
        fail(f"BH {pid}: invalid returned packet mass {dm}")
    if not (math.isfinite(q) and q > 0):
        fail(f"BH {pid}: invalid kick root q={q}")
    if np_ < MIN_TARGETS or nm_ < MIN_TARGETS:
        fail(f"BH {pid}: packet fired with too few targets ({np_},{nm_})")
    if not (math.isfinite(de) and de <= 3.1e-8):
        fail(f"BH {pid}: kinetic packet closure error {de}")
    if not (math.isfinite(pb) and pb <= 3.1e-8):
        fail(f"BH {pid}: bipolar kick momentum imbalance {pb}")

snaps = sorted(glob.glob(os.path.join(OUTPUT, "snap_*.hdf5")))
if not snaps:
    fail(f"no snapshots found in {OUTPUT}")

seen_bh = False
seen_threshold = False
seen_jet_buffer = False
hist = {}

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
        required = [
            "Masses", "BH_Mass", "BH_DiskMass", "BH_WindBufferMass",
            "BH_MdotFeed", "BH_MdotHorizon", "BH_MdotWind",
            "BH_EWind", "BH_EJet", "BH_EthWind", "BH_EthJet",
        ]
        missing = [name for name in required if name not in g]
        if missing:
            fail(f"{os.path.basename(path)} missing fields {missing}")

        ids = np.asarray(g["ParticleIDs"][:], dtype=np.uint64)
        for i, pid0 in enumerate(ids):
            pid = int(pid0)
            dyn = float(g["Masses"][i])
            bh = float(g["BH_Mass"][i])
            disk = float(g["BH_DiskMass"][i])
            wbuf = float(g["BH_WindBufferMass"][i])
            proc = float(g["BH_MdotFeed"][i])
            mdoth = float(g["BH_MdotHorizon"][i])
            mdotw = float(g["BH_MdotWind"][i])
            ewind = float(g["BH_EWind"][i])
            ejet = float(g["BH_EJet"][i])
            ethw = float(g["BH_EthWind"][i])
            ethj = float(g["BH_EthJet"][i])

            vals = [dyn, bh, disk, wbuf, proc, mdoth, mdotw, ewind, ejet, ethw, ethj]
            if not np.all(np.isfinite(vals)) or min(vals) < 0:
                fail(f"{os.path.basename(path)} BH {pid}: invalid feedback field")

            if not np.isclose(dyn, bh + disk + wbuf, rtol=6e-6, atol=1e-12):
                fail(f"{os.path.basename(path)} BH {pid}: dynamical mass ledger failed")
            if not np.isclose(proc, mdoth + mdotw, rtol=6e-6, atol=1e-12):
                fail(f"{os.path.basename(path)} BH {pid}: inner mass partition failed")

            if ethw > 0:
                seen_threshold = True
                if not np.isclose(ethj / ethw, JET_FACTOR / WIND_FACTOR, rtol=2e-5):
                    fail(f"{os.path.basename(path)} BH {pid}: wind/jet threshold ratio mismatch")

            if ejet > 0:
                seen_jet_buffer = True

            recs = hist.setdefault(pid, [])
            if recs and ejet + 1e-12 < recs[-1]["ejet"]:
                fail(f"{os.path.basename(path)} BH {pid}: jet buffer decreased during wind-only Iteration 7")
            recs.append({"ejet": ejet, "ewind": ewind, "wbuf": wbuf})

            print(
                f"  BH {pid}: Mwindbuf={wbuf*MASS_TO_MSUN:.6e} Msun "
                f"Ewind={ewind:.6e} Ejet={ejet:.6e} "
                f"EthWind={ethw:.6e} EthJet={ethj:.6e}"
            )
            print(
                f"    proc={proc*RATE_TO_MSUN_YR:.6e}, "
                f"H={mdoth*RATE_TO_MSUN_YR:.6e}, "
                f"wind={mdotw*RATE_TO_MSUN_YR:.6e} Msun/yr"
            )

if not seen_bh:
    fail("no BH formed")
if not seen_threshold:
    fail("no positive feedback binding threshold was recorded")
if not seen_jet_buffer:
    fail("jet energy never accumulated; wind-only/jet-buffer separation was not exercised")

print()
print("PASS: FFR-MACER Iteration-7 bipolar wind-feedback test")
print(f"  resolved wind packets fired : {len(events)}")
print("  every logged packet passed exact kinetic-energy closure")
print("  every logged packet passed zero-net-kick-momentum closure")
print("  both wind lobes satisfied the minimum target count")
print("  returned wind mass remains in the BH/gas dynamical mass ledger")
print("  jet energy remains buffered and is not consumed in Iteration 7")
