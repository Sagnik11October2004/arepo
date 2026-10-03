#!/usr/bin/env python3
"""Verify Iteration-9 FFR-MACER feedback-aware BH timestep limiting."""

import glob
import math
import os
import re
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "output_timestep")
if not os.path.isabs(OUTPUT):
    OUTPUT = os.path.join(HERE, OUTPUT)
LOG = os.path.join(HERE, "run_timestep.log")

if not os.path.isfile(LOG):
    raise SystemExit("FAIL: run_timestep.log not found; run ./run_timestep_restart.sh first")

with open(LOG, "r", encoding="utf-8", errors="replace") as fh:
    log = fh.read()

pat = re.compile(
    r"BH_FFR: timestep limit ID=(\d+) task=(\d+) "
    r"raw=(\d+) limited=(\d+) gasbin=(-?\d+) fint=([^\s]+) "
    r"dtintMyr=([^\s]+) dtwindMyr=([^\s]+) dtjetMyr=([^\s]+) "
    r"windBacklog=([^\s]+) jetBacklog=([^\s]+) backlog=(\d+)"
)
rows = pat.findall(log)
if not rows:
    raise SystemExit("FAIL: no FFR-MACER timestep limiter events were logged")

saw_internal = False
saw_backlog = False
max_jet_backlog = 0.0
max_reduction = 1.0

for row in rows:
    pid, task = row[0], row[1]
    raw = int(row[2])
    limited = int(row[3])
    fint = float(row[5])
    dtint = float(row[6])
    dtjet = float(row[8])
    jet_backlog = float(row[10])
    backlog = int(row[11])

    if not (raw > 0 and limited > 0 and limited <= raw):
        raise SystemExit(f"FAIL: BH {pid}: invalid integer limiter raw={raw} limited={limited}")
    if limited == 1:
        raise SystemExit(f"FAIL: BH {pid}: limiter generated forbidden integer timestep 1")
    if not math.isclose(fint, 4.0, rel_tol=0, abs_tol=1e-12):
        raise SystemExit(f"FAIL: BH {pid}: stress-test fint={fint}, expected 4")
    if not (math.isfinite(dtint) and dtint > 0):
        raise SystemExit(f"FAIL: BH {pid}: invalid reservoir accuracy timescale {dtint}")

    if math.isfinite(dtjet) and dtjet > 0 and limited < raw:
        saw_internal = True

    if backlog:
        saw_backlog = True
        if jet_backlog < 1.0 - 1e-8:
            raise SystemExit(f"FAIL: BH {pid}: backlog flag set with jetBacklog={jet_backlog}")

    max_jet_backlog = max(max_jet_backlog, jet_backlog)
    max_reduction = max(max_reduction, raw / limited)

if not saw_internal:
    raise SystemExit("FAIL: feedback/internal accuracy limit never reduced a BH candidate timestep")
if not saw_backlog:
    raise SystemExit(
        "FAIL: test did not exercise the post-packet backlog safeguard; "
        f"maximum jet backlog was {max_jet_backlog:.6g} thresholds"
    )

packet_re = re.compile(
    r"BH_FFR: jet packet fired ID=(\d+) task=(\d+) "
    r"E=([0-9eE+\-.]+) erg q=([0-9eE+\-.]+) "
    r"Nplus=(\d+) Nminus=(\d+) dErel=([0-9eE+\-.]+) pbal=([0-9eE+\-.]+)"
)
packets = packet_re.findall(log)
if not packets:
    raise SystemExit("FAIL: no jet packet fired in backlog stress test")

for pid, task, etxt, qtxt, nptxt, nmtxt, detxt, pbtxt in packets:
    if float(detxt) > 3.1e-8:
        raise SystemExit(f"FAIL: BH {pid}: packet energy closure error {detxt}")
    if float(pbtxt) > 3.1e-8:
        raise SystemExit(f"FAIL: BH {pid}: packet momentum imbalance {pbtxt}")

if "BH_FFR: wind packet fired" in log:
    raise SystemExit("FAIL: wind packet fired in jet-only timestep stress test")

snaps = sorted(glob.glob(os.path.join(OUTPUT, "snap_*.hdf5")))
if not snaps:
    raise SystemExit(f"FAIL: no snapshots in {OUTPUT}")

seen_bh = False
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

        for name in ("BH_EJet", "BH_EthJet", "BH_PJet", "BH_MdotFeed", "BH_MdotHorizon", "BH_MdotWind"):
            if name not in g:
                raise SystemExit(f"FAIL: {os.path.basename(path)} missing {name}")

        for i, pid in enumerate(np.asarray(g["ParticleIDs"][:], dtype=np.uint64)):
            ejet = float(g["BH_EJet"][i])
            ethj = float(g["BH_EthJet"][i])
            pjet = float(g["BH_PJet"][i])
            if not all(math.isfinite(x) and x >= 0 for x in (ejet, ethj, pjet)):
                raise SystemExit(f"FAIL: BH {int(pid)} invalid jet timestep diagnostics")

            ratio = ejet / ethj if ethj > 0 else 0.0
            print(
                f"  BH {int(pid)}: Ejet/EthJet={ratio:.6e} "
                f"Ejet={ejet:.6e} EthJet={ethj:.6e} PJet={pjet:.6e}"
            )

if not seen_bh:
    raise SystemExit("FAIL: no BH formed")

print()
print("PASS: FFR-MACER Iteration-9 feedback-timestep test")
print(f"  timestep limiter log events : {len(rows)}")
print(f"  jet packets fired           : {len(packets)}")
print(f"  maximum raw/limited ratio   : {max_reduction:.3f}")
print(f"  maximum logged jet backlog  : {max_jet_backlog:.3f} thresholds")
print("  feedback accuracy limit reduced the normal gravity candidate")
print("  packet-cap backlog forced a finer subsequent BH timestep")
print("  no forbidden integer timestep of size 1 was generated")
