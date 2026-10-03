#!/usr/bin/env python3
"""Verify the Iteration-6 inner accretion/energetics transaction."""

import glob
import math
import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "output_inner")
if not os.path.isabs(OUTPUT):
    OUTPUT = os.path.join(HERE, OUTPUT)

HUBBLE = 0.68
UNIT_MASS_MSUN = 1.0e10
UNIT_MASS_G = 1.989e43
UNIT_LENGTH_CM = 3.085678e21
UNIT_VELOCITY_CM_S = 1.0e5
UNIT_TIME_S = UNIT_LENGTH_CM / UNIT_VELOCITY_CM_S
SEC_PER_YEAR = 365.25 * 24.0 * 3600.0
C_CGS = 2.99792458e10
C_CODE = C_CGS / UNIT_VELOCITY_CM_S
EPS = 0.1
RHOT = 1200.0

MASS_TO_MSUN = UNIT_MASS_MSUN / HUBBLE
RATE_TO_MSUN_YR = UNIT_MASS_MSUN / (UNIT_TIME_S / SEC_PER_YEAR)
POWER_TO_CGS = UNIT_MASS_G * UNIT_VELOCITY_CM_S**2 / UNIT_TIME_S
ENERGY_TO_CGS = UNIT_MASS_G * UNIT_VELOCITY_CM_S**2 / HUBBLE

MODE = {0: "ADIOS", 1: "TRUNCATED", 2: "COLD"}


def fail(msg):
    raise SystemExit("FAIL: " + msg)


def close(a, b, rtol=3e-5, atol=1e-12):
    return np.allclose(a, b, rtol=rtol, atol=atol)


def cold_drive_luminosity_code(mdoth, mdotedd):
    if mdoth <= 0:
        return 0.0
    dotmh = mdoth / mdotedd
    if dotmh <= 1.0:
        return EPS * mdoth * C_CODE**2
    ledd = EPS * mdotedd * C_CODE**2
    return ledd * (1.0 + math.log(dotmh))


def radiative_luminosity_code(mdoth, mdotedd):
    if mdoth <= 0:
        return 0.0
    dotmh = mdoth / mdotedd
    if dotmh <= 1.0:
        return EPS * mdoth * C_CODE**2 * min(1.0, 10.0 * dotmh)
    ledd = EPS * mdotedd * C_CODE**2
    return ledd * (1.0 + math.log(dotmh))


snaps = sorted(glob.glob(os.path.join(OUTPUT, "snap_*.hdf5")))
if not snaps:
    fail(f"no snapshots found in {OUTPUT}")

seen_bh = False
seen_processed = False
seen_horizon = False
seen_wind = False
seen_hot_power = False
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
            "BH_MdotSupply", "BH_MdotFeed", "BH_MdotEdd",
            "BH_ProcessedEddRatio", "BH_MdotHorizon", "BH_MdotWind",
            "BH_Mode", "BH_EWind", "BH_EJet", "BH_Lbol", "BH_PWind",
            "BH_PJet",
        ]
        missing = [x for x in required if x not in g]
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
            edd = float(g["BH_MdotEdd"][i])
            ratio = float(g["BH_ProcessedEddRatio"][i])
            mdoth = float(g["BH_MdotHorizon"][i])
            mdotw = float(g["BH_MdotWind"][i])
            mode = int(g["BH_Mode"][i])
            ewind = float(g["BH_EWind"][i])
            ejet = float(g["BH_EJet"][i])
            lbol = float(g["BH_Lbol"][i])
            pwind = float(g["BH_PWind"][i])
            pjet = float(g["BH_PJet"][i])

            vals = [dyn, bh, disk, wbuf, proc, edd, ratio, mdoth, mdotw, ewind, ejet, lbol, pwind, pjet]
            if not np.all(np.isfinite(vals)):
                fail(f"{os.path.basename(path)} BH {pid}: non-finite inner-flow field")
            if min(vals) < 0:
                fail(f"{os.path.basename(path)} BH {pid}: negative inner-flow field")
            if mode not in MODE:
                fail(f"{os.path.basename(path)} BH {pid}: invalid mode {mode}")

            if not close(dyn, bh + disk + wbuf, rtol=5e-6):
                fail(f"{os.path.basename(path)} BH {pid}: dynamical mass ledger failed")
            if not close(proc, mdoth + mdotw, rtol=5e-6):
                fail(f"{os.path.basename(path)} BH {pid}: MdotFeed != MdotHorizon + MdotWind")
            if not (edd > 0):
                fail(f"{os.path.basename(path)} BH {pid}: invalid MdotEdd")
            if not close(ratio, proc / edd, rtol=5e-6, atol=1e-9):
                fail(f"{os.path.basename(path)} BH {pid}: processed Eddington ratio mismatch")

            expected_lbol = radiative_luminosity_code(mdoth, edd)
            if not close(lbol, expected_lbol, rtol=7e-5, atol=1e-12):
                fail(f"{os.path.basename(path)} BH {pid}: luminosity closure mismatch")

            if proc > 0:
                seen_processed = True
                if mdoth > 0:
                    seen_horizon = True
                if mdotw > 0:
                    seen_wind = True

                if mode in (0, 1):
                    q = 0.02 / ratio
                    rtr = 3.0 * q * q
                    r0 = max(3.0, min(rtr, RHOT))
                    fh = math.sqrt(3.0 / r0)
                    if not close(mdoth, fh * proc, rtol=7e-5):
                        fail(f"{os.path.basename(path)} BH {pid}: hot retention mismatch")

                    vwind = 0.2 * C_CODE / math.sqrt(2.0 * r0)
                    expected_pwind = 0.5 * mdotw * vwind * vwind
                    expected_pjet = 0.125 * mdoth * C_CODE * C_CODE
                    if not close(pwind, expected_pwind, rtol=8e-5):
                        fail(f"{os.path.basename(path)} BH {pid}: hot wind power mismatch")
                    if not close(pjet, expected_pjet, rtol=8e-5):
                        fail(f"{os.path.basename(path)} BH {pid}: hot jet power mismatch")
                    if pwind > 0 and pjet > 0:
                        seen_hot_power = True
                else:
                    ldrive_code = cold_drive_luminosity_code(mdoth, edd)
                    ldrive_cgs = ldrive_code * POWER_TO_CGS
                    expected_w_msun_yr = 0.28 * (ldrive_cgs / 1.0e45) ** 0.85
                    expected_w = expected_w_msun_yr / RATE_TO_MSUN_YR
                    if not close(mdotw, expected_w, rtol=1e-4, atol=1e-12):
                        fail(f"{os.path.basename(path)} BH {pid}: cold Gofford mass-loss closure mismatch")

                    v_kms = min(2.5e4 * (ldrive_cgs / 1.0e45) ** 0.4, 1.0e5)
                    v_code = v_kms * 1.0e5 / UNIT_VELOCITY_CM_S
                    expected_pwind = 0.5 * mdotw * v_code * v_code
                    if not close(pwind, expected_pwind, rtol=1e-4):
                        fail(f"{os.path.basename(path)} BH {pid}: cold wind power mismatch")
                    if abs(pjet) > 1e-12:
                        fail(f"{os.path.basename(path)} BH {pid}: cold state has non-zero jet power")

            recs = hist.setdefault(pid, [])
            if recs:
                prev = recs[-1]
                if wbuf + 1e-12 < prev["wbuf"]:
                    fail(f"{os.path.basename(path)} BH {pid}: wind mass buffer decreased before injection stage")
                if ewind + 1e-12 < prev["ewind"] or ejet + 1e-12 < prev["ejet"]:
                    fail(f"{os.path.basename(path)} BH {pid}: mechanical energy buffer decreased before injection stage")
                if bh + 1e-12 < prev["bh"]:
                    fail(f"{os.path.basename(path)} BH {pid}: BH mass decreased")

            recs.append({"bh": bh, "wbuf": wbuf, "ewind": ewind, "ejet": ejet})

            print(
                f"  BH {pid}: mode={MODE[mode]} "
                f"Mdisk={disk*MASS_TO_MSUN:.6e} Msun "
                f"Mwindbuf={wbuf*MASS_TO_MSUN:.6e} Msun"
            )
            print(
                f"    proc={proc*RATE_TO_MSUN_YR:.6e}, "
                f"H={mdoth*RATE_TO_MSUN_YR:.6e}, "
                f"wind={mdotw*RATE_TO_MSUN_YR:.6e} Msun/yr, "
                f"dotm_proc={ratio:.6e}"
            )
            print(
                f"    Lbol={lbol*POWER_TO_CGS:.6e} erg/s "
                f"Pwind={pwind*POWER_TO_CGS:.6e} erg/s "
                f"Pjet={pjet*POWER_TO_CGS:.6e} erg/s"
            )
            print(
                f"    Ewind={ewind*ENERGY_TO_CGS:.6e} erg "
                f"Ejet={ejet*ENERGY_TO_CGS:.6e} erg"
            )

if not seen_bh:
    fail("no BH formed")
if not seen_processed:
    fail("no non-zero reservoir processing was observed")
if not seen_horizon:
    fail("no horizon flow was observed")
if not seen_wind:
    fail("no wind flow was observed")
if not seen_hot_power:
    fail("no hot-state wind+jet power was observed")

print()
print("PASS: FFR-MACER Iteration-6 inner-flow test")
print("  exact processed = horizon + wind rate partition")
print("  dynamical mass ledger is consistent")
print("  hot retention/wind/jet formulae match")
print("  luminosity closure matches")
print("  wind/jet buffers accumulate without injection")
print("  no hard Eddington mass cap is assumed by the verifier")
