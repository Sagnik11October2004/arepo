#!/usr/bin/env python3
from __future__ import annotations

import math
import sys
from pathlib import Path

import h5py
import numpy as np

OMEGA_M = 0.31
OMEGA_B = 0.048
H = 0.68
L_MPC_H = 0.5
N = 128 ** 3
ZSTART = 99.0


def main() -> None:
    path = Path(sys.argv[1] if len(sys.argv) > 1 else "ics.hdf5")
    with h5py.File(path, "r") as f:
        h = f["Header"].attrs
        time = float(h["Time"])
        redshift = float(h["Redshift"])
        box = float(h["BoxSize"])
        counts = np.asarray(h["NumPart_Total"], dtype=np.uint64)
        mass_table = np.asarray(h["MassTable"], dtype=np.float64)

        if counts[1] != N:
            raise SystemExit(f"Expected {N} parent Type-1 particles; found {counts[1]}")
        if counts.sum() != N:
            raise SystemExit(f"Expected a DM-only MUSIC IC with {N} particles; counts={counts.tolist()}")
        if not math.isclose(redshift, ZSTART, rel_tol=0.0, abs_tol=2.0e-6):
            raise SystemExit(f"Expected z={ZSTART}; Header/Redshift={redshift}")
        if not math.isclose(time, 1.0 / (1.0 + ZSTART), rel_tol=0.0, abs_tol=2.0e-8):
            raise SystemExit(f"Unexpected Header/Time={time}")

    # rho_crit,0 = 2.77536627e11 h^2 Msun/Mpc^3.
    rho_crit = 2.77536627e11 * H * H
    volume_mpc3 = (L_MPC_H / H) ** 3
    mtot = rho_crit * OMEGA_M * volume_mpc3
    mgas = mtot * OMEGA_B / OMEGA_M / N
    mdm = mtot * (OMEGA_M - OMEGA_B) / OMEGA_M / N

    print(f"IC OK: {path}")
    print(f"  z_start          = {redshift:.6f}")
    print(f"  Header BoxSize   = {box:.9g} code length")
    print(f"  parent particles = {N:,}")
    print(f"  MassTable        = {mass_table.tolist()}")
    print(f"  expected post-split m_gas = {mgas:.3f} Msun")
    print(f"  expected post-split m_DM  = {mdm:.3f} Msun")
    print(f"  expected total box mass   = {mtot:.6e} Msun")
    if not mgas < 1.0e4:
        raise SystemExit("Requested gas mass resolution <1e4 Msun is not satisfied")


if __name__ == "__main__":
    main()
