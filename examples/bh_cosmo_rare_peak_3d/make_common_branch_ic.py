#!/usr/bin/env python3
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import h5py
import numpy as np

H = 0.68
EXPECTED_Z = 20.0
EXPECTED_BH_MASS_MSUN = 1.0e5


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: make_common_branch_ic.py SOURCE_SNAPSHOT DEST_IC")
    src = Path(sys.argv[1]).resolve()
    dst = Path(sys.argv[2]).resolve()
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)

    with h5py.File(dst, "r+") as f:
        z = float(f["Header"].attrs["Redshift"])
        if abs(z - EXPECTED_Z) > 2.0e-3:
            raise SystemExit(f"Expected common branch state at z~20, found z={z}")

        if "PartType5" not in f:
            raise SystemExit("Common z=20 state has no PartType5 BH")
        g = f["PartType5"]
        nbh = len(g["ParticleIDs"])
        if nbh != 1:
            raise SystemExit(f"Expected exactly one BH, found {nbh}")

        m_code = float(np.asarray(g["Masses"])[0])
        m_msun = m_code * 1.0e10 / H
        if abs(m_msun / EXPECTED_BH_MASS_MSUN - 1.0) > 2.0e-6:
            raise SystemExit(
                f"Expected untouched {EXPECTED_BH_MASS_MSUN:.6g} Msun seed; found {m_msun:.9g} Msun"
            )

        # Snapshot BH_* fields are output-only in the current FFR I/O table.
        # Remove them explicitly so every branch is unmistakably a fresh,
        # zero-reservoir Type-5 IC initialized from the same phase-space state.
        for name in list(g.keys()):
            if name.startswith("BH_"):
                del g[name]

    print(f"Wrote common fresh branch IC: {dst}")
    print(f"  z={z:.6f}, N_BH=1, M_BH={m_msun:.6g} Msun")


if __name__ == "__main__":
    main()
