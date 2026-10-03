#!/usr/bin/env python3
"""Inspect FFR-MACER BH evolution in the FoF seeding smoke-test snapshots."""

import glob
import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = sys.argv[1] if len(sys.argv) > 1 else "output"
if not os.path.isabs(OUTPUT):
    OUTPUT = os.path.join(HERE, OUTPUT)
SNAPS = sorted(glob.glob(os.path.join(OUTPUT, "snap_*.hdf5")))

HUBBLE = 0.68
UNIT_MASS_MSUN = 1.0e10
UNIT_LENGTH_CM = 3.085678e21
UNIT_VELOCITY_CM_S = 1.0e5
SEC_PER_YEAR = 365.25 * 24.0 * 3600.0
UNIT_TIME_S = UNIT_LENGTH_CM / UNIT_VELOCITY_CM_S

MASS_TO_MSUN = UNIT_MASS_MSUN / HUBBLE
RATE_TO_MSUN_YR = UNIT_MASS_MSUN / (UNIT_TIME_S / SEC_PER_YEAR)

FIELDS = [
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
    "BH_Coherence",
    "BH_EWind",
    "BH_EJet",
    "BH_Lbol",
    "BH_PWind",
    "BH_PJet",
    "BH_SigmaDM",
    "BH_EthWind",
    "BH_EthJet",
]


def scalar(group, name, i):
    if name not in group:
        return None
    x = np.asarray(group[name][i])
    return float(x) if x.ndim == 0 else x


if not SNAPS:
    sys.exit(f"No snapshots found in {OUTPUT}.")

history = {}

print("FFR-MACER black-hole snapshot diagnostic")
print("=" * 72)

for path in SNAPS:
    with h5py.File(path, "r") as f:
        z = float(f["Header"].attrs["Redshift"])
        a = float(f["Header"].attrs["Time"])
        g = f.get("PartType5")
        nbh = 0 if g is None else len(g["ParticleIDs"])

        print(f"\n{os.path.basename(path)}: a={a:.9g}, z={z:.6f}, N_BH={nbh}")
        if not nbh:
            continue

        ids = np.asarray(g["ParticleIDs"][:], dtype=np.uint64)
        for i, pid_raw in enumerate(ids):
            pid = int(pid_raw)
            bh = scalar(g, "BH_Mass", i)
            dyn = scalar(g, "Masses", i)
            disk = scalar(g, "BH_DiskMass", i)
            wbuf = scalar(g, "BH_WindBufferMass", i)

            print(f"  BH ID {pid}")
            if bh is not None:
                print(f"    physical BH mass : {bh * MASS_TO_MSUN:.9e} Msun")
            if dyn is not None:
                print(f"    dynamical mass   : {dyn * MASS_TO_MSUN:.9e} Msun")
            if disk is not None:
                print(f"    reservoir mass   : {disk * MASS_TO_MSUN:.9e} Msun")
            if wbuf is not None:
                print(f"    wind-buffer mass : {wbuf * MASS_TO_MSUN:.9e} Msun")

            for name in ("BH_MdotSupply", "BH_MdotFeed", "BH_MdotEdd", "BH_MdotHorizon", "BH_MdotWind"):
                val = scalar(g, name, i)
                if val is not None:
                    suffix = " (processed reservoir rate)" if name == "BH_MdotFeed" else ""
                    print(f"    {name:16s}: {val * RATE_TO_MSUN_YR:.9e} Msun/yr{suffix}")

            for name in ("BH_ProcessedEddRatio", "BH_Mode", "BH_Coherence", "BH_EWind", "BH_EJet",
                         "BH_Lbol", "BH_PWind", "BH_PJet", "BH_SigmaDM", "BH_EthWind", "BH_EthJet"):
                val = scalar(g, name, i)
                if val is not None:
                    print(f"    {name:16s}: {val}")

            for name in ("BH_DiscAxis", "BH_JetAxis"):
                val = scalar(g, name, i)
                if val is not None:
                    print(f"    {name:16s}: {np.asarray(val)}")

            history.setdefault(pid, []).append(
                {
                    "snapshot": os.path.basename(path),
                    "z": z,
                    "mass": bh,
                    "dyn": dyn,
                    "disk": disk,
                    "wbuf": wbuf,
                }
            )

print("\n" + "=" * 72)
print("Evolution summary")

for pid, recs in sorted(history.items()):
    first, last = recs[0], recs[-1]
    dm = (last["mass"] - first["mass"]) * MASS_TO_MSUN
    ddyn = (last["dyn"] - first["dyn"]) * MASS_TO_MSUN
    ddisk = (last["disk"] - first["disk"]) * MASS_TO_MSUN
    dwbuf = (last["wbuf"] - first["wbuf"]) * MASS_TO_MSUN

    print(f"  BH ID {pid}: {first['snapshot']} -> {last['snapshot']}")
    print(f"    Delta M_BH        = {dm:.9e} Msun")
    print(f"    Delta M_dyn       = {ddyn:.9e} Msun")
    print(f"    Delta M_reservoir = {ddisk:.9e} Msun")
    print(f"    Delta M_windbuf   = {dwbuf:.9e} Msun")

    if abs(dm) <= 0.1 and abs(ddisk) <= 0.1 and abs(dwbuf) <= 0.1:
        print("    RESULT            = no measurable inner mass transfer")
    else:
        print("    RESULT            = inner-flow mass transfer detected; verify with verify_inner.py")
