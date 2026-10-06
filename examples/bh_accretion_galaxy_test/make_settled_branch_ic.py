#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent


def copy_required(src: h5py.File, dst: h5py.File, ptype: int) -> int:
    name = f"PartType{ptype}"
    if name not in src:
        return 0
    sg = src[name]
    dg = dst.create_group(name)
    for key in ("Coordinates", "Velocities", "ParticleIDs", "Masses"):
        dg.create_dataset(key, data=np.asarray(sg[key]))
    if ptype == 0:
        dg.create_dataset("InternalEnergy", data=np.asarray(sg["InternalEnergy"]))
    return len(sg["ParticleIDs"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("snapshot", type=Path)
    ap.add_argument("--output", type=Path, default=HERE / "ics" / "settled_galaxy.hdf5")
    args = ap.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.snapshot, "r") as src, h5py.File(args.output, "w") as dst:
        counts = np.zeros(6, dtype=np.uint32)
        for ptype in (0, 1, 4, 5):
            counts[ptype] = copy_required(src, dst, ptype)

        h = dst.create_group("Header")
        h.attrs["NumPart_ThisFile"] = counts
        h.attrs["NumPart_Total"] = counts
        h.attrs["NumPart_Total_HighWord"] = np.zeros(6, dtype=np.uint32)
        h.attrs["MassTable"] = np.zeros(6, dtype=np.float64)
        h.attrs["Time"] = 0.0
        h.attrs["Redshift"] = 0.0
        h.attrs["BoxSize"] = float(src["Header"].attrs["BoxSize"])
        h.attrs["NumFilesPerSnapshot"] = 1
        h.attrs["Omega0"] = 0.0
        h.attrs["OmegaLambda"] = 0.0
        h.attrs["HubbleParam"] = 1.0
        h.attrs["Flag_Sfr"] = 0
        h.attrs["Flag_Cooling"] = 0
        h.attrs["Flag_StellarAge"] = 0
        h.attrs["Flag_Metals"] = 0
        h.attrs["Flag_Feedback"] = 0
        h.attrs["Flag_DoublePrecision"] = 1
        h.attrs["Flag_Entropy_ICs"] = 0

        meta = dst.create_group("ICMetadata")
        meta.attrs["Generator"] = "bh_accretion_galaxy_test/make_settled_branch_ic.py"
        meta.attrs["SourceSnapshot"] = str(args.snapshot)
        meta.attrs["ClockReset"] = 1
        meta.attrs["BHSubgridStateReset"] = 1

    print(f"Wrote fresh post-settling IC: {args.output}")
    print(f"Counts: gas={counts[0]}, DM={counts[1]}, stars={counts[4]}, BH={counts[5]}")


if __name__ == "__main__":
    main()
