#!/usr/bin/env python3
"""Create a deterministic two-BH snapshot for the Iteration-11 merger test."""

import json
import os
import shutil
import sys

import h5py
import numpy as np

if len(sys.argv) != 3:
    raise SystemExit("usage: make_merger_snapshot.py SOURCE DEST")

src, dst = sys.argv[1:3]
if not os.path.isfile(src):
    raise SystemExit(f"missing source snapshot: {src}")

shutil.copy2(src, dst)

with h5py.File(dst, "r+") as f:
    if "PartType5" not in f or len(f["PartType5"]["ParticleIDs"]) != 1:
        raise SystemExit("source snapshot must contain exactly one Type-5 BH")

    a = float(f["Header"].attrs["Time"])
    box = float(f["Header"].attrs["BoxSize"])
    g = f["PartType5"]

    max_id = 0
    for key in f:
        if key.startswith("PartType") and "ParticleIDs" in f[key]:
            ids = np.asarray(f[key]["ParticleIDs"][:], dtype=np.uint64)
            if len(ids):
                max_id = max(max_id, int(ids.max()))

    original = {}
    attrs = {}
    for name, ds in list(g.items()):
        arr = ds[...]
        if arr.shape[0] != 1:
            raise SystemExit(f"unexpected PartType5 dataset shape {name}: {arr.shape}")
        original[name] = arr
        attrs[name] = dict(ds.attrs)

    for name in list(g.keys()):
        del g[name]

    for name, arr in original.items():
        doubled = np.concatenate([arr, arr[:1]], axis=0)
        ds = g.create_dataset(name, data=doubled, dtype=arr.dtype)
        for key, val in attrs[name].items():
            ds.attrs[key] = val

    id0 = int(g["ParticleIDs"][0])
    id1 = max_id + 1
    g["ParticleIDs"][1] = np.array(id1, dtype=g["ParticleIDs"].dtype)

    pos0 = np.asarray(g["Coordinates"][0], dtype=float)
    vel0 = np.asarray(g["Velocities"][0], dtype=float)
    m0 = float(g["Masses"][0])

    # The canonical smoke test has R_acc=0.10 proper kpc/h. Put the duplicate
    # only 0.05 proper units away, comfortably inside the 2 R_acc overlap rule.
    proper_sep = 0.05
    coord_sep = proper_sep / a
    pos1 = pos0.copy()
    pos1[0] = (pos1[0] + coord_sep) % box
    g["Coordinates"][1] = pos1.astype(g["Coordinates"].dtype)

    vel1 = vel0.copy()
    vel1[0] += 10.0
    g["Velocities"][1] = vel1.astype(g["Velocities"].dtype)

    # Keep equal dynamical masses so the exact merger velocity is the midpoint.
    g["Masses"][1] = np.array(m0, dtype=g["Masses"].dtype)

    # Keep output-only BH mass visually consistent in the synthetic file if
    # present. RestartFlag=2 still reconstructs FFR state from Masses.
    if "BH_Mass" in g:
        g["BH_Mass"][1] = np.array(m0, dtype=g["BH_Mass"].dtype)
    for name in ("BH_DiskMass", "BH_WindBufferMass", "BH_EWind", "BH_EJet"):
        if name in g:
            g[name][1] = np.zeros_like(g[name][1])

    nthis = np.asarray(f["Header"].attrs["NumPart_ThisFile"]).copy()
    ntot = np.asarray(f["Header"].attrs["NumPart_Total"]).copy()
    nthis[5] += 1
    ntot[5] += 1
    f["Header"].attrs.modify("NumPart_ThisFile", nthis)
    f["Header"].attrs.modify("NumPart_Total", ntot)

expected = {
    "source": os.path.abspath(src),
    "a": a,
    "id0": id0,
    "id1": id1,
    "survivor_id": min(id0, id1),
    "mass0": m0,
    "mass1": m0,
    "merged_mass": 2.0 * m0,
    "vel0": vel0.tolist(),
    "vel1": vel1.tolist(),
    "merged_vel": (0.5 * (vel0 + vel1)).tolist(),
    "proper_separation": proper_sep,
    "merge_threshold": 0.20,
}

sidecar = os.path.splitext(dst)[0] + "_expected.json"
with open(sidecar, "w", encoding="utf-8") as fh:
    json.dump(expected, fh, indent=2)

print(f"Created {dst}")
print(json.dumps(expected, indent=2))
