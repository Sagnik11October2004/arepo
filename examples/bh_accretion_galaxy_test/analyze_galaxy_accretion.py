#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from pathlib import Path
import re

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
LOG = HERE / "logs" / "galaxy_env_shell.log"
OUT = HERE / "outputs" / "galaxy_env_shell"
RATE_CODE_TO_MSUN_YR = 10.2271202634
MYR_PER_CODE = 977.792354298

RAW_KEYS = ("TNG", "BOOSTED", "AM", "FFR_VOLUME", "FFR_SHELL", "FFR_ENV")
ENV_KEYS = ("FENV", "FTH", "FWIND", "FROT", "thermalX", "windX", "rotX",
            "velcoh", "Mcoh", "Rsh", "cssh", "vbulksh", "vphish", "vdynsh")

def extract(line: str, key: str) -> float | None:
    m = re.search(rf'(?<![A-Za-z0-9_]){re.escape(key)}='
                  r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)', line)
    return float(m.group(1)) if m else None

def parse_log() -> list[dict]:
    if not LOG.exists():
        raise FileNotFoundError(LOG)
    rows = []
    with LOG.open("r", errors="replace") as f:
        for line in f:
            if "BH_BENCHMARK_ALL:" not in line or "FFR_ENV=" not in line:
                continue
            row = {}
            keys = ("time", "timeMyr", "captureon", "Mgas", "rho", "cs", "vrel",
                    "Vphi", "nH", "boost", "amlim", *RAW_KEYS, *ENV_KEYS)
            for key in keys:
                x = extract(line, key)
                if x is not None:
                    row[key] = x
            if "timeMyr" not in row:
                continue
            for key in RAW_KEYS:
                if key in row:
                    row[f"{key}_Msun_yr"] = row[key] * RATE_CODE_TO_MSUN_YR
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

def parse_snapshots() -> list[dict]:
    rows = []
    for path in sorted(OUT.glob("snap_*.hdf5")):
        with h5py.File(path, "r") as f:
            time_code = float(f["Header"].attrs["Time"])
            row = {"snapshot": path.name, "time_code": time_code,
                   "time_Myr": time_code * MYR_PER_CODE}
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

def median_in(rows, key, t0, t1):
    vals = [r[key] for r in rows if t0 <= r.get("timeMyr", -1) < t1 and key in r]
    return float(np.median(vals)) if vals else math.nan

def print_rate_summary(rows):
    if not rows:
        print("No FFR_ENV diagnostic rows found.")
        return
    windows = [("settled/pre-on", 15.0, 20.0),
               ("early-on", 20.0, 30.0),
               ("late", 50.0, 80.0001)]
    print()
    print("MEDIAN RAW ACCRETION ESTIMATORS [Msun/yr]")
    print("=" * 112)
    print(f"{'window':18s}" + "".join(f"{k:>15s}" for k in RAW_KEYS) + f"{'FENV':>12s}")
    print("-" * 112)
    for name, t0, t1 in windows:
        vals = [median_in(rows, f"{k}_Msun_yr", t0, t1) for k in RAW_KEYS]
        fenv = median_in(rows, "FENV", t0, t1)
        print(f"{name:18s}" + "".join(f"{v:15.4e}" for v in vals) + f"{fenv:12.4f}")
    on = [r for r in rows if r.get("captureon", 0) > 0.5]
    off = [r for r in rows if r.get("captureon", 0) < 0.5]
    if off and on:
        print()
        print(f"BH-off diagnostics: {len(off)} samples; BH-on diagnostics: {len(on)} samples")
        print(f"First capture-on diagnostic: {min(r['timeMyr'] for r in on):.6f} Myr")

def make_plots(rows, states):
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"Plotting skipped: {e}")
        return
    if rows:
        fig, ax = plt.subplots()
        t = np.array([r["timeMyr"] for r in rows])
        for key in RAW_KEYS:
            y = np.array([r.get(f"{key}_Msun_yr", np.nan) for r in rows])
            ax.plot(t, y, label=key)
        ax.axvline(20.0, linestyle="--")
        ax.set_yscale("log")
        ax.set_xlabel("Time [Myr]")
        ax.set_ylabel("Raw accretion estimator [Msun/yr]")
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        fig.savefig(HERE / "galaxy_accretion_rates.png", dpi=180)
        plt.close(fig)

        fig, ax = plt.subplots()
        for key in ("FENV", "FTH", "FWIND", "FROT"):
            y = np.array([r.get(key, np.nan) for r in rows])
            ax.plot(t, y, label=key)
        ax.axvline(20.0, linestyle="--")
        ax.set_xlabel("Time [Myr]")
        ax.set_ylabel("Environmental factor")
        ax.set_ylim(0.0, 1.05)
        ax.legend()
        fig.tight_layout()
        fig.savefig(HERE / "galaxy_environment_factors.png", dpi=180)
        plt.close(fig)

    if states:
        fig, ax = plt.subplots()
        t = np.array([r["time_Myr"] for r in states])
        for key in ("BH_MdotSupply_Msun_yr", "BH_MdotFeed_Msun_yr",
                    "BH_MdotHorizon_Msun_yr", "BH_MdotWind_Msun_yr"):
            y = np.array([r.get(key, np.nan) for r in states])
            ax.plot(t, y, marker="o", label=key.replace("_Msun_yr", ""))
        ax.axvline(20.0, linestyle="--")
        ax.set_yscale("symlog", linthresh=1.0e-12)
        ax.set_xlabel("Time [Myr]")
        ax.set_ylabel("Selected-path rate [Msun/yr]")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(HERE / "galaxy_bh_reservoir_rates.png", dpi=180)
        plt.close(fig)

def main():
    rows = parse_log()
    states = parse_snapshots()
    write_csv(HERE / "galaxy_accretion_estimators.csv", rows)
    write_csv(HERE / "galaxy_bh_state.csv", states)
    print_rate_summary(rows)
    print()
    print(f"Raw estimator samples: {len(rows)}")
    print(f"BH snapshot states: {len(states)}")
    print("Saved galaxy_accretion_estimators.csv")
    print("Saved galaxy_bh_state.csv")
    make_plots(rows, states)
    print("Saved comparison plots when matplotlib is available.")

if __name__ == "__main__":
    main()
