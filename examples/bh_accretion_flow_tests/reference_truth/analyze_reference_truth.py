#!/usr/bin/env python3
"""Measure and convergence-test the high-resolution reference accretion rates."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
FLOW = HERE.parent

UNIT_MASS_MSUN = 1.0e10
MYR_PER_CODE = 977.792354298
YR_PER_CODE = MYR_PER_CODE * 1.0e6
KM_CGS = 1.0e5
PC_CGS = 3.0856775814913673e18
YR_CGS = 365.25 * 86400.0

LEVEL_ORDER = {"HR4": 4, "HR5": 5, "HR6": 6}
FLUX_RADII_PC = (0.5, 1.0, 2.0, 4.0, 8.0, 12.5)


def scalar_dataset(group, name: str) -> float:
    a = np.asarray(group[name])
    return float(a.reshape(-1)[0])


def snapshot_files(output_dir: Path) -> list[Path]:
    return sorted(output_dir.glob("snap_*.hdf5"))


def read_bh_history(output_dir: Path) -> dict:
    snaps = snapshot_files(output_dir)
    if len(snaps) < 4:
        raise RuntimeError(f"{output_dir}: need >=4 snapshots, found {len(snaps)}")

    time_yr = []
    cumulative_msun = []
    instant_msun_yr = []

    for snap in snaps:
        with h5py.File(snap, "r") as f:
            tcode = float(f["Header"].attrs["Time"])
            bh = f["PartType5"]
            time_yr.append(tcode * YR_PER_CODE)
            cumulative_msun.append(
                scalar_dataset(bh, "BH_BenchmarkCumulativeBHMassGrowth") * UNIT_MASS_MSUN
            )
            if "BH_BenchmarkMdotRealized" in bh:
                instant_msun_yr.append(
                    scalar_dataset(bh, "BH_BenchmarkMdotRealized")
                    * UNIT_MASS_MSUN / YR_PER_CODE
                )
            else:
                instant_msun_yr.append(float("nan"))

    t = np.asarray(time_yr)
    m = np.asarray(cumulative_msun)
    inst = np.asarray(instant_msun_yr)

    order = np.argsort(t)
    t = t[order]
    m = m[order]
    inst = inst[order]

    # Use the final 40% of the evolution to suppress initial transients.
    late = t >= (0.60 * t[-1])
    if np.count_nonzero(late) < 4:
        late = np.arange(t.size) >= max(0, t.size - 4)

    tl = t[late]
    ml = m[late]

    if tl.size < 2 or tl[-1] <= tl[0]:
        raise RuntimeError(f"{output_dir}: insufficient late-time history")

    # Linear slope of cumulative swallowed mass is the primary mean rate.
    slope = float(np.polyfit(tl, ml, 1)[0])

    dt = np.diff(tl)
    dm = np.diff(ml)
    interval_rates = dm / dt
    interval_mean = float(np.mean(interval_rates))
    interval_std = float(np.std(interval_rates, ddof=1)) if interval_rates.size > 1 else 0.0
    variability = interval_std / max(abs(interval_mean), 1.0e-300)

    inst_late = inst[late]
    finite_inst = inst_late[np.isfinite(inst_late)]
    instant_mean = float(np.mean(finite_inst)) if finite_inst.size else float("nan")

    return {
        "snapshot_count": len(snaps),
        "time_final_yr": float(t[-1]),
        "late_start_yr": float(tl[0]),
        "cumulative_growth_final_msun": float(m[-1]),
        "mean_sink_mdot_msun_yr": slope,
        "interval_mean_mdot_msun_yr": interval_mean,
        "interval_std_mdot_msun_yr": interval_std,
        "late_variability_fraction": variability,
        "instantaneous_realized_mean_msun_yr": instant_mean,
    }


def periodic_delta(delta: np.ndarray, box_kpc: float) -> np.ndarray:
    return delta - box_kpc * np.rint(delta / box_kpc)


def shell_flux_from_snapshot(snap: Path, radius_pc: float, finest_dx_pc: float) -> float:
    # Do not pretend to measure a shell whose radius itself is unresolved.
    if radius_pc < 4.0 * finest_dx_pc:
        return float("nan")

    with h5py.File(snap, "r") as f:
        box_kpc = float(f["Header"].attrs["BoxSize"])
        gas = f["PartType0"]
        bh = f["PartType5"]

        x = np.asarray(gas["Coordinates"], dtype=np.float64)
        v = np.asarray(gas["Velocities"], dtype=np.float64)
        m_msun = np.asarray(gas["Masses"], dtype=np.float64) * UNIT_MASS_MSUN

        xb = np.asarray(bh["Coordinates"], dtype=np.float64).reshape(-1, 3)[0]
        vb = np.asarray(bh["Velocities"], dtype=np.float64).reshape(-1, 3)[0]

    dr_kpc = periodic_delta(x - xb[None, :], box_kpc)
    dr_pc = 1000.0 * dr_kpc
    r_pc = np.linalg.norm(dr_pc, axis=1)

    # Use a shell wide enough to contain several finest cells while remaining
    # local on larger radii.
    full_width_pc = max(4.0 * finest_dx_pc, 0.20 * radius_pc)
    half = 0.5 * full_width_pc
    mask = (r_pc >= radius_pc - half) & (r_pc < radius_pc + half)
    if np.count_nonzero(mask) < 8:
        return float("nan")

    rhat = dr_pc[mask] / r_pc[mask, None]
    vr_kms = np.sum((v[mask] - vb[None, :]) * rhat, axis=1)

    conv = KM_CGS / PC_CGS * YR_CGS
    # Positive means net inward.
    mdot = -np.sum(m_msun[mask] * vr_kms) / full_width_pc * conv
    return float(mdot)


def late_fluxes(output_dir: Path, finest_dx_pc: float) -> dict:
    snaps = snapshot_files(output_dir)
    if not snaps:
        return {}

    use = snaps[-min(5, len(snaps)):]
    result = {}
    for radius in FLUX_RADII_PC:
        vals = [
            shell_flux_from_snapshot(s, radius, finest_dx_pc)
            for s in use
        ]
        vals = np.asarray([x for x in vals if math.isfinite(x)])
        key = str(radius).replace(".", "p")
        result[f"flux_r{key}pc_mean_msun_yr"] = (
            float(np.mean(vals)) if vals.size else float("nan")
        )
        result[f"flux_r{key}pc_std_msun_yr"] = (
            float(np.std(vals, ddof=1)) if vals.size > 1 else
            (0.0 if vals.size == 1 else float("nan"))
        )
    return result


def metadata_for(run: dict) -> dict:
    p = (
        HERE / run["family"] / "ics" / run["level"]
        / f"truth_{run['family']}_{run['tag']}_{run['level']}.metadata.json"
    )
    return json.loads(p.read_text())


def analyze_completed_runs(manifest: dict) -> list[dict]:
    rows = []

    for run in manifest["runs"]:
        output_dir = FLOW / run["output_dir"]
        if not (output_dir / "end").exists():
            continue

        try:
            hist = read_bh_history(output_dir)
        except Exception as exc:
            print(f"WARNING: {run['family']} {run['tag']} {run['level']}: {exc}")
            continue

        meta = metadata_for(run)
        flux = late_fluxes(output_dir, float(run["finest_dx_pc"]))

        row = {
            "family": run["family"],
            "tag": run["tag"],
            "control_value": run["control_value"],
            "level": run["level"],
            "finest_dx_pc": run["finest_dx_pc"],
            "sink_radius_pc": run["sink_radius_pc"],
            "absorber_A": run["absorber_A"],
            "absorber_max_fraction": run["absorber_max_fraction"],
            **hist,
            **flux,
        }

        if run["family"] == "bhl":
            row["analytic_bhl_mdot_msun_yr"] = meta["analytic_bhl_mdot_msun_per_yr"]
            row["r_bhl_pc"] = meta["r_bhl_pc"]
        else:
            row["nominal_supply_mdot_msun_yr"] = meta["mdot_true_msun_per_yr"]

        rows.append(row)

    rows.sort(key=lambda r: (r["family"], r["control_value"], LEVEL_ORDER[r["level"]]))
    return rows


def build_case_summary(rows: list[dict]) -> list[dict]:
    grouped: dict[tuple, dict[str, dict]] = {}
    for row in rows:
        key = (row["family"], row["tag"], row["control_value"])
        grouped.setdefault(key, {})[row["level"]] = row

    summary = []
    for (family, tag, value), levels in sorted(grouped.items()):
        item = {
            "family": family,
            "tag": tag,
            "control_value": value,
            "completed_levels": sorted(levels, key=lambda x: LEVEL_ORDER[x]),
        }

        for level in ("HR4", "HR5", "HR6"):
            if level in levels:
                item[f"{level}_mdot_msun_yr"] = levels[level]["mean_sink_mdot_msun_yr"]
                item[f"{level}_variability_fraction"] = levels[level]["late_variability_fraction"]
                item[f"{level}_sink_radius_pc"] = levels[level]["sink_radius_pc"]

        if "HR5" in levels and "HR6" in levels:
            a = levels["HR5"]["mean_sink_mdot_msun_yr"]
            b = levels["HR6"]["mean_sink_mdot_msun_yr"]
            conv = abs(b - a) / max(abs(b), 1.0e-300)
            item["HR5_to_HR6_relative_change"] = conv
            item["converged_10pct"] = bool(conv <= 0.10)
            item["candidate_HR6_actual_mdot_msun_yr"] = b
            item["actual_mdot_msun_yr"] = b if conv <= 0.10 else None
        else:
            item["HR5_to_HR6_relative_change"] = None
            item["converged_10pct"] = False
            item["candidate_HR6_actual_mdot_msun_yr"] = (
                levels["HR6"]["mean_sink_mdot_msun_yr"] if "HR6" in levels else None
            )
            item["actual_mdot_msun_yr"] = None

        if family == "bhl" and "HR6" in levels:
            analytic = levels["HR6"].get("analytic_bhl_mdot_msun_yr")
            item["analytic_bhl_mdot_msun_yr"] = analytic
            if analytic and item["candidate_HR6_actual_mdot_msun_yr"] is not None:
                item["HR6_over_analytic_BHL"] = (
                    item["candidate_HR6_actual_mdot_msun_yr"] / analytic
                )

        summary.append(item)

    return summary


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                fields.append(key)
                seen.add(key)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    manifest_path = HERE / "reference_run_manifest.json"
    if not manifest_path.exists():
        raise SystemExit(
            "Missing reference_run_manifest.json. Run "
            "python3 reference_truth/generate_reference_runs.py first."
        )

    manifest = json.loads(manifest_path.read_text())
    rows = analyze_completed_runs(manifest)
    summary = build_case_summary(rows)

    (HERE / "reference_truth_summary.json").write_text(
        json.dumps({"runs": rows, "cases": summary}, indent=2, sort_keys=True) + "\n"
    )
    write_csv(HERE / "reference_truth_runs.csv", rows)
    write_csv(HERE / "reference_truth_cases.csv", summary)

    print()
    print("HIGH-RESOLUTION REFERENCE ACCRETION RATES")
    print("=" * 108)
    print(
        f"{'CASE':20s} {'HR4':>12s} {'HR5':>12s} {'HR6':>12s} "
        f"{'d56':>10s} {'status':>12s}"
    )
    print("-" * 108)

    for s in summary:
        case = f"{s['family']}:{s['tag']}"
        def fmt(level):
            x = s.get(f"{level}_mdot_msun_yr")
            return f"{x:.4e}" if x is not None else "-"
        d = s.get("HR5_to_HR6_relative_change")
        dstr = f"{d:.3%}" if d is not None else "-"
        status = "CONVERGED" if s.get("converged_10pct") else "NOT YET"
        print(
            f"{case:20s} {fmt('HR4'):>12s} {fmt('HR5'):>12s} {fmt('HR6'):>12s} "
            f"{dstr:>10s} {status:>12s}"
        )

    print()
    print(f"Completed reference runs analyzed: {len(rows)} / {len(manifest['runs'])}")
    print("A case is promoted to actual_mdot only when HR5->HR6 changes by <=10%.")
    print("For turbulent cases, late-time variability is reported separately rather than used as a rejection criterion.")
    print("Saved reference_truth/reference_truth_summary.json")


if __name__ == "__main__":
    main()
