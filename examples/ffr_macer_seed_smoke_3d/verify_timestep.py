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

HUBBLE = 0.68
UNIT_MASS_MSUN = 1.0e10
UNIT_LENGTH_CM = 3.085678e21
UNIT_VELOCITY_CM_S = 1.0e5
SEC_PER_YEAR = 365.25 * 24.0 * 3600.0
UNIT_TIME_S = UNIT_LENGTH_CM / UNIT_VELOCITY_CM_S
RATE_TO_MSUN_YR = UNIT_MASS_MSUN / (UNIT_TIME_S / SEC_PER_YEAR)

if not os.path.isfile(LOG):
    raise SystemExit("FAIL: run_timestep.log not found; run ./run_timestep_restart.sh first")

with open(LOG, "r", encoding="utf-8", errors="replace") as fh:
    log = fh.read()

pat = re.compile(
    r"BH_FFR: timestep limit ID=(\d+) task=(\d+) "
    r"raw=(\d+) limited=(\d+) gasbin=(-?\d+) fint=([^\s]+) "
    r"dtintMyr=([^\s]+) dtwindMyr=([^\s]+) dtjetMyr=([^\s]+) "
    r"accuracyLimited=(\d+) backlogLimit=(\d+) "
    r"windBacklog=([^\s]+) jetBacklog=([^\s]+) backlog=(\d+) backlogApplied=(\d+)"
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
    accuracy_limited = int(row[9])
    backlog_limit = int(row[10])
    jet_backlog = float(row[12])
    backlog = int(row[13])
    backlog_applied = int(row[14])

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

    if not (2 <= accuracy_limited <= raw):
        raise SystemExit(
            f"FAIL: BH {pid}: invalid accuracy-limited step {accuracy_limited} from raw {raw}"
        )

    if backlog_applied:
        saw_backlog = True
        if not backlog:
            raise SystemExit(f"FAIL: BH {pid}: backlogApplied without backlog flag")
        if jet_backlog < 1.0 - 1e-8:
            raise SystemExit(f"FAIL: BH {pid}: backlog limiter applied with jetBacklog={jet_backlog}")
        expected_backlog = max(2, accuracy_limited // 2)
        if backlog_limit != expected_backlog or limited != expected_backlog:
            raise SystemExit(
                f"FAIL: BH {pid}: backlog limiter recursively/incorrectly reduced step "
                f"accuracy={accuracy_limited} backlogLimit={backlog_limit} limited={limited} "
                f"expected={expected_backlog}"
            )

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

        for name in ("BH_EJet", "BH_EthJet", "BH_PJet", "BH_MdotSupply", "BH_MdotFeed", "BH_MdotHorizon", "BH_MdotWind"):
            if name not in g:
                raise SystemExit(f"FAIL: {os.path.basename(path)} missing {name}")

        for i, pid in enumerate(np.asarray(g["ParticleIDs"][:], dtype=np.uint64)):
            ejet = float(g["BH_EJet"][i])
            ethj = float(g["BH_EthJet"][i])
            pjet = float(g["BH_PJet"][i])
            supply = float(g["BH_MdotSupply"][i])
            if not all(math.isfinite(x) and x >= 0 for x in (ejet, ethj, pjet, supply)):
                raise SystemExit(f"FAIL: BH {int(pid)} invalid jet timestep diagnostics")

            supply_phys = supply * RATE_TO_MSUN_YR
            if supply_phys > 1.0e-2:
                raise SystemExit(
                    f"FAIL: BH {int(pid)} spurious supply-rate spike {supply_phys:.6e} Msun/yr; "
                    "capture rate must use gas hydro timesteps, not the shorter BH feedback timestep"
                )

            ratio = ejet / ethj if ethj > 0 else 0.0
            print(
                f"  BH {int(pid)}: Ejet/EthJet={ratio:.6e} "
                f"Ejet={ejet:.6e} EthJet={ethj:.6e} PJet={pjet:.6e} "
                f"MdotSupply={supply_phys:.6e} Msun/yr"
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
print("  packet-cap backlog forced exactly one bin below the accuracy-limited candidate")
print("  backlog limiting is anchored and does not recursively collapse the BH timebin")
print("  resolved-supply diagnostics use contributing gas hydro timesteps")
print("  no forbidden integer timestep of size 1 was generated")
