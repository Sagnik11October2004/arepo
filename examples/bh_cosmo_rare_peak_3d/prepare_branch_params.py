#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import h5py

HERE = Path(__file__).resolve().parent
BASE = HERE / "param_seed_window.txt"
COMMON_IC = HERE / "ics" / "common_z20_seeded.hdf5"

MODELS = {
    0: "tng_bondi",
    1: "boosted_bondi",
    2: "am_bondi",
    3: "ffr_volume",
    4: "ffr_shell",
    5: "convj_shell_ffr",
}
FEEDBACKS = {
    "none": (0, None),
    "tng": (1, 1),      # TNG feedback requires direct target.
    "macer": (2, 0),    # MACER feedback requires reservoir target.
}


def set_param(text: str, key: str, value: str) -> str:
    lines = text.splitlines()
    found = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped and not stripped.startswith("%") and stripped.split()[0] == key:
            lines[i] = f"{key:<40s}{value}"
            found = True
            break
    if not found:
        raise KeyError(key)
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--z-end", type=float, default=15.0)
    ap.add_argument("--none-target", choices=("direct", "reservoir"), default="direct")
    args = ap.parse_args()

    if not (0.0 < args.z_end < 20.0):
        raise SystemExit("--z-end must satisfy 0 < z_end < 20")
    if not COMMON_IC.exists():
        raise SystemExit(f"Missing common branch IC: {COMMON_IC}. Run ./run_to_z20.sh first.")

    with h5py.File(COMMON_IC, "r") as f:
        a0 = float(f["Header"].attrs["Time"])
        z0 = float(f["Header"].attrs["Redshift"])
        nbh = len(f["PartType5"]["ParticleIDs"]) if "PartType5" in f else 0
    if abs(z0 - 20.0) > 2.0e-3 or nbh != 1:
        raise SystemExit(f"Expected one-BH common IC at z=20; found z={z0}, N_BH={nbh}")

    amax = 1.0 / (1.0 + args.z_end)
    base = BASE.read_text()
    params_dir = HERE / "params"
    params_dir.mkdir(exist_ok=True)

    # Common snapshot schedule for the default z=20 -> z_end branches.
    candidate_z = [19, 18, 17, 16, 15, 14, 12, 10, 8, 6]
    with (HERE / "output_branches.txt").open("w") as f:
        for z in candidate_z:
            if args.z_end < z < 20:
                f.write(f"{1.0/(1.0+z):.17g} 2\n")
        f.write(f"{amax:.17g} 2\n")

    none_target = 1 if args.none_target == "direct" else 0

    for model, mname in MODELS.items():
        for fbname, (fbmodel, forced_target) in FEEDBACKS.items():
            target = none_target if forced_target is None else forced_target
            text = base
            updates = {
                "InitCondFile": "./ics/common_z20_seeded",
                "OutputDir": f"./outputs/{mname}__{fbname}",
                "OutputListFilename": "./output_branches.txt",
                "TimeBegin": f"{a0:.17g}",
                "TimeMax": f"{amax:.17g}",
                "MaxSizeTimestep": "0.001",
                "BHSeedMinRedshift": "100.0",
                "BHFreeFallA": "1.0e-3",
                "BHBenchmarkAccretionModel": str(model),
                "BHBenchmarkAccretionTarget": str(target),
                "BHBenchmarkFeedbackModel": str(fbmodel),
                "BHBenchmarkEddingtonFactor": "1.0",
            }
            for key, value in updates.items():
                text = set_param(text, key, value)
            (params_dir / f"{mname}__{fbname}.txt").write_text(text)

    print(f"Prepared {len(MODELS)*len(FEEDBACKS)} branch parameter files in {params_dir}")
    print(f"Common start: z={z0:.6f}; branch end: z={args.z_end:g}")
    print(f"No-feedback target: {args.none_target}")
    print("Target rules: TNG feedback=direct; MACER feedback=reservoir.")
    print()
    print("IMPORTANT: mode 5 is available, but its current benchmark implementation")
    print("uses the unit-efficiency geometric shell base. Freeze/calibrate its final")
    print("A_ff/C_ff normalization before treating the mode-5 cosmological branch as")
    print("a production science result.")


if __name__ == "__main__":
    main()
