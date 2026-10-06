#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from pathlib import Path
import re

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
RATE_CODE_TO_MSUN_YR = 10.2271202634
MYR_PER_CODE = 977.792354298
SETTLE_MYR = 20.0

MODELS = {
    0: ("tng_bondi", "TNG"),
    1: ("boosted_bondi", "BOOSTED"),
    2: ("am_bondi", "AM"),
    3: ("ffr_volume", "FFR_VOLUME"),
    4: ("ffr_shell", "FFR_SHELL"),
    5: ("ffr_env", "FFR_ENV"),
}
RAW_KEYS = ("TNG", "BOOSTED", "AM", "FFR_VOLUME", "FFR_SHELL", "FFR_ENV")
ENV_KEYS = ("FENV", "FTH", "FWIND", "FROT", "thermalX", "windX", "rotX",
            "velcoh", "Mcoh", "Rsh", "cssh", "vbulksh", "vphish", "vdynsh")


def extract(line: str, key: str) -> float | None:
    m = re.search(
        rf'(?<![A-Za-z0-9_]){re.escape(key)}='
        r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)',
        line,
    )
    return float(m.group(1)) if m else None


def parse_log(path: Path, stage: str, model_id: int | None, model_name: str | None) -> list[dict]:
    if not path.exists():
        return []
    rows = []
    with path.open("r", errors="replace") as f:
        for line in f:
            if "BH_BENCHMARK_ALL:" not in line or "FFR_ENV=" not in line:
                continue
            row = {"stage": stage, "model_id": model_id, "model_name": model_name}
            keys = ("time", "timeMyr", "captureon", "Mgas", "rho", "cs", "vrel",
                    "Vphi", "nH", "boost", "amlim", *RAW_KEYS, *ENV_KEYS)
            for key in keys:
                x = extract(line, key)
                if x is not None:
                    row[key] = x
            if "timeMyr" not in row:
                continue
            row["galaxyTimeMyr"] = row["timeMyr"] if stage == "settle" else SETTLE_MYR + row["timeMyr"]
            for key in RAW_KEYS:
                if key in row:
                    row[f"{key}_Msun_yr"] = row[key] * RATE_CODE_TO_MSUN_YR
            if model_id is not None:
                selected_key = MODELS[model_id][1]
                if selected_key in row:
                    row["selected_raw_Msun_yr"] = row[selected_key] * RATE_CODE_TO_MSUN_YR
            rows.append(row)
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields, seen = [], set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fields.append(key)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def flatten_value(row: dict, key: str, value):
    arr = np.asarray(value)
    if arr.ndim == 0:
        row[key] = float(arr)
    elif arr.size == 1:
        row[key] = float(arr.reshape(-1)[0])
    else:
        for i, x in enumerate(arr.reshape(-1)):
            row[f"{key}_{i}"] = float(x)


def parse_snapshot_dir(path: Path, stage: str, model_id: int | None, model_name: str | None) -> list[dict]:
    rows = []
    for snap in sorted(path.glob("snap_*.hdf5")):
        with h5py.File(snap, "r") as f:
            time_code = float(f["Header"].attrs["Time"])
            time_myr = time_code * MYR_PER_CODE
            row = {
                "stage": stage,
                "model_id": model_id,
                "model_name": model_name,
                "snapshot": snap.name,
                "time_code": time_code,
                "time_Myr": time_myr,
                "galaxyTimeMyr": time_myr if stage == "settle" else SETTLE_MYR + time_myr,
            }
            if "PartType5" not in f:
                continue
            bh = f["PartType5"]
            for key in sorted(bh.keys()):
                if key.startswith("BH_") or key in ("Coordinates", "Velocities", "Masses"):
                    flatten_value(row, key, bh[key][0])
            for key in ("BH_Mass", "BH_DiskMass", "BH_WindBufferMass", "Masses"):
                if key in row:
                    row[f"{key}_Msun"] = row[key] * 1.0e10
            for key in ("BH_MdotSupply", "BH_MdotFeed", "BH_MdotEdd",
                        "BH_MdotHorizon", "BH_MdotWind",
                        "BH_BenchmarkMdotRaw", "BH_BenchmarkMdotOperational",
                        "BH_BenchmarkMdotRealized"):
                if key in row:
                    row[f"{key}_Msun_yr"] = row[key] * RATE_CODE_TO_MSUN_YR
            rows.append(row)
    return rows


def median(vals):
    vals = [v for v in vals if np.isfinite(v)]
    return float(np.median(vals)) if vals else math.nan


def summarize(raw_rows, state_rows):
    settle = [r for r in raw_rows if r["stage"] == "settle" and 15 <= r["galaxyTimeMyr"] <= 20]
    print()
    print("COMMON SETTLED ENVIRONMENT: MEDIAN RAW ESTIMATORS AT 15-20 MYR [Msun/yr]")
    print("=" * 98)
    for key in RAW_KEYS:
        vals = [r.get(f"{key}_Msun_yr", math.nan) for r in settle]
        print(f"{key:12s} {median(vals):.6e}")

    print()
    print("POST-SETTLE BRANCH COMPARISON")
    print("=" * 124)
    print(f"{'model':18s} {'selected raw':>14s} {'realized':>14s} {'Mdot_H':>14s} {'Mdisk_final':>14s} {'MBH_final':>14s}")
    print("-" * 124)
    for model_id, (name, selected_key) in MODELS.items():
        rr = [r for r in raw_rows if r.get("model_id") == model_id and 50 <= r["galaxyTimeMyr"] <= 80]
        sr = [r for r in state_rows if r.get("model_id") == model_id]
        selected = median([r.get("selected_raw_Msun_yr", math.nan) for r in rr])
        realized = median([r.get("BH_BenchmarkMdotRealized_Msun_yr", math.nan) for r in sr if 50 <= r["galaxyTimeMyr"] <= 80])
        horizon = median([r.get("BH_MdotHorizon_Msun_yr", math.nan) for r in sr if 50 <= r["galaxyTimeMyr"] <= 80])
        final = max(sr, key=lambda r: r["galaxyTimeMyr"]) if sr else {}
        mdisk = final.get("BH_DiskMass_Msun", math.nan)
        mbh = final.get("BH_Mass_Msun", final.get("Masses_Msun", math.nan))
        print(f"{name:18s} {selected:14.4e} {realized:14.4e} {horizon:14.4e} {mdisk:14.4e} {mbh:14.4e}")


def make_plots(raw_rows, state_rows):
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"Plotting skipped: {e}")
        return

    fig, ax = plt.subplots()
    for model_id, (name, _) in MODELS.items():
        rr = [r for r in raw_rows if r.get("model_id") == model_id]
        if not rr:
            continue
        ax.plot([r["galaxyTimeMyr"] for r in rr],
                [r.get("selected_raw_Msun_yr", np.nan) for r in rr],
                label=name)
    ax.axvline(SETTLE_MYR, linestyle="--")
    ax.set_yscale("log")
    ax.set_xlabel("Galaxy time [Myr]")
    ax.set_ylabel("Selected raw supply estimator [Msun/yr]")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(HERE / "galaxy_branch_selected_rates.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots()
    for model_id, (name, _) in MODELS.items():
        sr = [r for r in state_rows if r.get("model_id") == model_id]
        if not sr:
            continue
        ax.plot([r["galaxyTimeMyr"] for r in sr],
                [r.get("BH_Mass_Msun", r.get("Masses_Msun", np.nan)) for r in sr],
                marker="o", label=name)
    ax.set_xlabel("Galaxy time [Myr]")
    ax.set_ylabel("BH mass [Msun]")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(HERE / "galaxy_branch_bh_mass.png", dpi=180)
    plt.close(fig)

    env = [r for r in raw_rows if r.get("model_id") == 5]
    if env:
        fig, ax = plt.subplots()
        for key in ("FENV", "FTH", "FWIND", "FROT"):
            ax.plot([r["galaxyTimeMyr"] for r in env], [r.get(key, np.nan) for r in env], label=key)
        ax.set_xlabel("Galaxy time [Myr]")
        ax.set_ylabel("Mode-5 environmental factor")
        ax.set_ylim(0.0, 1.05)
        ax.legend()
        fig.tight_layout()
        fig.savefig(HERE / "galaxy_environment_factors.png", dpi=180)
        plt.close(fig)


def main():
    raw_rows = parse_log(HERE / "logs" / "settle.log", "settle", None, None)
    state_rows = parse_snapshot_dir(HERE / "outputs" / "settle", "settle", None, None)

    for model_id, (name, _) in MODELS.items():
        raw_rows += parse_log(HERE / "logs" / f"model_{model_id}_{name}.log", "branch", model_id, name)
        state_rows += parse_snapshot_dir(HERE / "outputs" / f"model_{model_id}_{name}", "branch", model_id, name)

    write_csv(HERE / "galaxy_accretion_estimators.csv", raw_rows)
    write_csv(HERE / "galaxy_bh_state.csv", state_rows)
    summarize(raw_rows, state_rows)
    make_plots(raw_rows, state_rows)
    print()
    print(f"Saved {len(raw_rows)} raw-estimator samples and {len(state_rows)} BH snapshot states.")
    print("Outputs: galaxy_accretion_estimators.csv, galaxy_bh_state.csv, galaxy_branch_selected_rates.png, galaxy_branch_bh_mass.png")


if __name__ == "__main__":
    main()
