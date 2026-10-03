#!/usr/bin/env python3
import glob
import math
import os
import re
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = os.path.join(HERE, "output")
LOG = os.path.join(HERE, "run.log")

EXPECTED_PER_COMPONENT = 128 ** 3
HUBBLE = 0.68
UNIT_MASS_MSUN = 1.0e10
EXPECTED_SEED_MSUN = 1.0e5
ZMIN = 20.0


def fail(msg):
    print("FAIL:", msg)
    sys.exit(1)


def particle_ids(f, ptype):
    gname = f"PartType{ptype}"
    if gname not in f:
        return np.empty(0, dtype=np.uint64)
    return np.asarray(f[gname]["ParticleIDs"][:], dtype=np.uint64)


def component_mass(f, ptype):
    gname = f"PartType{ptype}"
    if gname not in f:
        return 0.0
    group = f[gname]
    n = len(group["ParticleIDs"])
    if "Masses" in group:
        return float(np.sum(group["Masses"][:], dtype=np.float64))
    table = np.asarray(f["Header"].attrs["MassTable"], dtype=np.float64)
    return float(n * table[ptype])


if not os.path.isfile(LOG):
    fail("run.log not found; run ./run.sh first")

with open(LOG, "r", encoding="utf-8", errors="replace") as fh:
    log = fh.read()

seed_re = re.compile(
    r"BH_FFR: seeded FoF group\s+(\d+)\s+at z=([0-9eE+\-.]+).*?MBH=([0-9eE+\-.]+)\s+Msun.*?ID=(\d+)"
)
events = seed_re.findall(log)
if not events:
    fail("no FoF seed event was reported in run.log")

for gr, ztxt, mtxt, pid in events:
    z = float(ztxt)
    m = float(mtxt)
    if not z > ZMIN:
        fail(f"seed group {gr} formed at z={z}, not strictly above z={ZMIN}")
    if not math.isclose(m, EXPECTED_SEED_MSUN, rel_tol=1.0e-8, abs_tol=1.0e-6):
        fail(f"seed group {gr} reported MBH={m} Msun instead of {EXPECTED_SEED_MSUN} Msun")

snaps = sorted(glob.glob(os.path.join(OUTPUT, "snap_*.hdf5")))
if len(snaps) < 2:
    fail("fewer than two HDF5 snapshots found")

records = []
reference_total_mass = None
first_bh_index = None
first_bh_ids = None

for isnap, path in enumerate(snaps):
    with h5py.File(path, "r") as f:
        z = float(f["Header"].attrs["Redshift"])
        ids0 = particle_ids(f, 0)
        ids1 = particle_ids(f, 1)
        ids5 = particle_ids(f, 5)

        n0, n1, n5 = len(ids0), len(ids1), len(ids5)
        if n0 + n5 != EXPECTED_PER_COMPONENT:
            fail(f"{os.path.basename(path)}: N_gas + N_BH = {n0+n5}, expected {EXPECTED_PER_COMPONENT}")
        if n1 != EXPECTED_PER_COMPONENT:
            fail(f"{os.path.basename(path)}: N_DM = {n1}, expected {EXPECTED_PER_COMPONENT}")

        total_mass = sum(component_mass(f, t) for t in range(6))
        if reference_total_mass is None:
            reference_total_mass = total_mass
        elif not math.isclose(total_mass, reference_total_mass, rel_tol=2.0e-6, abs_tol=1.0e-12):
            fail(
                f"{os.path.basename(path)}: snapshot mass changed from "
                f"{reference_total_mass:.12e} to {total_mass:.12e} code-mass units"
            )

        if n5 > 0:
            if first_bh_index is None:
                first_bh_index = isnap
                first_bh_ids = set(int(x) for x in ids5)

            g = f["PartType5"]
            required = [
                "Masses",
                "BH_Mass",
                "BH_DiskMass",
                "BH_WindBufferMass",
                "BH_MdotSupply",
                "BH_MdotFeed",
                "BH_MdotEdd",
                "BH_ProcessedEddRatio",
                "BH_MdotHorizon",
                "BH_MdotWind",
                "BH_Mode",
                "BH_Lbol",
                "BH_PWind",
                "BH_PJet",
            ]
            missing = [name for name in required if name not in g]
            if missing:
                fail(f"{os.path.basename(path)}: missing Type-5 fields {missing}")

            bh_mass_code = np.asarray(g["BH_Mass"][:], dtype=np.float64)
            bh_mass_msun = bh_mass_code * UNIT_MASS_MSUN / HUBBLE
            if not np.all(np.isfinite(bh_mass_msun)):
                fail(f"{os.path.basename(path)}: non-finite BH_Mass")
            if not np.allclose(bh_mass_msun, EXPECTED_SEED_MSUN, rtol=2.0e-6, atol=0.1):
                fail(f"{os.path.basename(path)}: BH_Mass is not the requested {EXPECTED_SEED_MSUN:g} Msun")

            dyn_mass = np.asarray(g["Masses"][:], dtype=np.float64)
            if not np.allclose(dyn_mass, bh_mass_code, rtol=2.0e-6, atol=1.0e-12):
                fail(f"{os.path.basename(path)}: seed dynamical mass differs from BH_Mass with capture disabled")

            disk = np.asarray(g["BH_DiskMass"][:], dtype=np.float64)
            windbuf = np.asarray(g["BH_WindBufferMass"][:], dtype=np.float64)
            mdot = np.asarray(g["BH_MdotSupply"][:], dtype=np.float64)
            mdot_feed = np.asarray(g["BH_MdotFeed"][:], dtype=np.float64)
            mdot_edd = np.asarray(g["BH_MdotEdd"][:], dtype=np.float64)
            edd_ratio = np.asarray(g["BH_ProcessedEddRatio"][:], dtype=np.float64)
            mode = np.asarray(g["BH_Mode"][:], dtype=np.int64)

            if np.any(np.abs(disk) > 1.0e-12):
                fail(f"{os.path.basename(path)}: non-zero reservoir mass although BHFreeFallA=0")
            if np.any(np.abs(windbuf) > 1.0e-12):
                fail(f"{os.path.basename(path)}: non-zero wind buffer before Iteration 6")
            if np.any(np.abs(mdot) > 1.0e-12):
                fail(f"{os.path.basename(path)}: non-zero supply rate although BHFreeFallA=0")
            if np.any(np.abs(mdot_feed) > 1.0e-12):
                fail(f"{os.path.basename(path)}: non-zero candidate feed rate for an empty reservoir")
            if np.any(~np.isfinite(mdot_edd)) or np.any(mdot_edd <= 0.0):
                fail(f"{os.path.basename(path)}: invalid BH_MdotEdd for seeded BH")
            if np.any(np.abs(edd_ratio) > 1.0e-12):
                fail(f"{os.path.basename(path)}: non-zero processed Eddington ratio for an empty reservoir")
            if np.any(mode != 0):
                fail(f"{os.path.basename(path)}: empty-reservoir BH should classify as ADIOS (mode 0)")

            for name in ("BH_MdotHorizon", "BH_MdotWind", "BH_Lbol", "BH_PWind", "BH_PJet"):
                arr = np.asarray(g[name][:], dtype=np.float64)
                if np.any(np.abs(arr) > 1.0e-12):
                    fail(f"{os.path.basename(path)}: {name} became non-zero before Iteration 6+")

        records.append((os.path.basename(path), z, n0, n1, n5, total_mass))

if first_bh_index is None:
    fail("run.log reported seeding but no snapshot contains a Type-5 BH")
if first_bh_index >= len(snaps) - 1:
    fail("the first BH appears only in the final snapshot; no post-seeding continuation snapshot exists")

with h5py.File(snaps[-1], "r") as f:
    final_ids = set(int(x) for x in particle_ids(f, 5))
if not first_bh_ids.issubset(final_ids):
    fail("one or more BH IDs from the first seeded snapshot are absent from the final snapshot")

print("PASS: FFR-MACER FoF seeding smoke test")
print(f"  seed events in log : {len(events)}")
print(f"  snapshots checked  : {len(records)}")
print(f"  first BH snapshot  : {records[first_bh_index][0]} (z={records[first_bh_index][1]:.6g})")
print(f"  final BH count     : {records[-1][4]}")
print(f"  conserved mass     : {reference_total_mass:.12e} code-mass units")
print("  gas/BH conversion  : N_gas + N_BH stayed exactly 128^3")
print("  DM count           : stayed exactly 128^3")
print("  seed mass          : 1.0e5 Msun")
print("  reservoir/supply   : zero as required for BHFreeFallA=0")
print("  Iteration-5 state  : empty reservoir -> zero candidate feed, ADIOS mode")
print("  later-stage physics: horizon/wind/radiation/feedback remain zero")
print("  post-seed survival : seeded BH IDs persist to the final snapshot")
