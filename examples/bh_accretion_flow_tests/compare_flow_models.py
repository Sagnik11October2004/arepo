#!/usr/bin/env python3
"""Compare all coarse accretion estimators on the same flow state.

Primary comparison is deliberately model-to-model.  Optional external
reference columns can be added with --with-reference, but they are not used
to normalize or calibrate mode 5.

Current mode 5:
    Mdot_5 = Mdot_shell,geom * f_M * f_j
    Cconv = max[(I-O)/(I+O), 0]
    M_eff = (1-Cconv) * M_bulk
    f_M = (1 + M_eff^2)^(-3/2)
    chi_j = j_shell^2 / (G M_cen R_acc)
    f_j = (1 + chi_j^2)^(-1/2)

I and O, shell bulk velocity, sound speed, and shell specific angular
momentum are all weighted by each cell's unit-efficiency shell-FFR
contribution m_i/t_ff,i.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
RATE_CODE_TO_MSUN_YR = 10.2271202634

FAMILIES = {
    "bhl": [
        ("M0p0", 0.0), ("M0p5", 0.5), ("M1p0", 1.0),
        ("M2p0", 2.0), ("M5p0", 5.0),
    ],
    "rotating": [
        ("f0p00", 0.0), ("f0p25", 0.25), ("f0p50", 0.50),
        ("f0p75", 0.75), ("f1p00", 1.00),
    ],
    "turbulent": [
        ("M0p5", 0.5), ("M1p0", 1.0), ("M2p0", 2.0), ("M5p0", 5.0),
    ],
}


def stem(family: str, tag: str, level: str) -> str:
    if family == "bhl":
        return f"bhl_{tag}_{level}"
    if family == "rotating":
        return f"rot_{tag}_{level}"
    if family == "turbulent":
        return f"turb_{tag}_{level}"
    raise ValueError(family)


def extract(line: str, key: str) -> float | None:
    m = re.search(
        rf'(?<![A-Za-z0-9_]){re.escape(key)}='
        r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)',
        line,
    )
    return float(m.group(1)) if m else None


def benchmark_lines(log: Path) -> list[str]:
    if not log.exists():
        return []
    out = []
    with log.open("r", errors="replace") as f:
        for line in f:
            if "BH_BENCHMARK_ALL:" in line and "FFR_CONVJ=" in line:
                out.append(line)
    return out


def read_aff(param: Path) -> float:
    if not param.exists():
        return 1.0e-3
    for line in param.read_text(errors="replace").splitlines():
        s = line.strip()
        if not s or s.startswith("%") or s.startswith("#"):
            continue
        fields = s.split()
        if len(fields) >= 2 and fields[0] == "BHFreeFallA":
            return float(fields[1])
    return 1.0e-3


def metadata(family: str, tag: str, level: str) -> dict:
    path = HERE / family / "ics" / level / f"{stem(family, tag, level)}.metadata.json"
    return json.loads(path.read_text())


def reference_rate(family: str, tag: str, value: float, level: str, tng_mdot: float) -> tuple[float, str]:
    if family == "bhl":
        if value < 1.0:
            return 0.25 * tng_mdot, "Bondi-gamma5over3"
        m = metadata(family, tag, level)
        return float(m["analytic_bhl_mdot_msun_per_yr"]), "BHL-baseline"
    m = metadata(family, tag, level)
    return float(m["mdot_true_msun_per_yr"]), "analytic-IC-supply"


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
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=("diagnostic", "evolved"), default="diagnostic")
    ap.add_argument("--sample", choices=("first", "last"), default="last")
    ap.add_argument("--with-reference", action="store_true")
    args = ap.parse_args()

    suffix = "diag" if args.stage == "diagnostic" else "evolve"
    rows: list[dict] = []
    missing: list[str] = []

    for family, cases in FAMILIES.items():
        for tag, value in cases:
            for level in ("L0", "L1", "L2"):
                case = stem(family, tag, level)
                log = HERE / family / "logs" / suffix / f"{case}.log"
                lines = benchmark_lines(log)
                if not lines:
                    missing.append(case)
                    continue
                line = lines[0] if args.sample == "first" else lines[-1]

                values = {}
                for key in (
                    "TNG", "BOOSTED", "AM", "FFR_VOLUME", "FFR_SHELL",
                    "FFR_SHELL_GEOM", "FFR_CONVJ", "FCORR", "FMACH", "FJ",
                    "CCONV", "Mbulk", "Meff", "chiJ", "rcirc",
                ):
                    x = extract(line, key)
                    if x is None:
                        raise ValueError(f"{case}: missing {key}=... in selected diagnostic line")
                    values[key] = x

                param = HERE / family / "params" / f"{case}_{suffix}.txt"
                aff = read_aff(param)
                geom = values["FFR_SHELL_GEOM"]
                if not (geom > 0):
                    continue

                row = {
                    "case": case,
                    "family": family,
                    "tag": tag,
                    "control_value": value,
                    "level": level,
                    "stage": args.stage,
                    "sample": args.sample,
                    "samples_available": len(lines),
                    "A_ff": aff,
                    "shell_geom_mdot_msun_yr": geom * RATE_CODE_TO_MSUN_YR,
                    "TNG_mdot_msun_yr": values["TNG"] * RATE_CODE_TO_MSUN_YR,
                    "BOOSTED_mdot_msun_yr": values["BOOSTED"] * RATE_CODE_TO_MSUN_YR,
                    "AM_mdot_msun_yr": values["AM"] * RATE_CODE_TO_MSUN_YR,
                    "FFR_VOLUME_mdot_msun_yr": values["FFR_VOLUME"] * RATE_CODE_TO_MSUN_YR,
                    "FFR_SHELL_mdot_msun_yr": values["FFR_SHELL"] * RATE_CODE_TO_MSUN_YR,
                    "FFR_CONVJ_mdot_msun_yr": values["FFR_CONVJ"] * RATE_CODE_TO_MSUN_YR,
                    "TNG_over_geom": values["TNG"] / geom,
                    "BOOSTED_over_geom": values["BOOSTED"] / geom,
                    "AM_over_geom": values["AM"] / geom,
                    "FFR_VOLUME_over_Aff_geom":
                        values["FFR_VOLUME"] / (aff * geom) if aff > 0 else math.nan,
                    "FFR_SHELL_over_Aff_geom":
                        values["FFR_SHELL"] / (aff * geom) if aff > 0 else math.nan,
                    "CONVJ_over_geom": values["FFR_CONVJ"] / geom,
                    "FCORR": values["FCORR"],
                    "FMACH": values["FMACH"],
                    "FJ": values["FJ"],
                    "CCONV": values["CCONV"],
                    "Mbulk": values["Mbulk"],
                    "Meff": values["Meff"],
                    "chiJ": values["chiJ"],
                    "rcirc_code_length": values["rcirc"],
                }

                if args.with_reference:
                    ref, kind = reference_rate(
                        family, tag, value, level, row["TNG_mdot_msun_yr"]
                    )
                    row["reference_kind"] = kind
                    row["reference_mdot_msun_yr"] = ref
                    for mode in ("TNG", "BOOSTED", "AM", "FFR_CONVJ"):
                        row[f"{mode}_over_reference"] = row[f"{mode}_mdot_msun_yr"] / ref
                    row["shell_geom_over_reference"] = row["shell_geom_mdot_msun_yr"] / ref

                rows.append(row)

    outfile = HERE / f"flow_model_comparison_{args.stage}_{args.sample}.csv"
    write_csv(outfile, rows)

    print()
    print("MODEL-TO-MODEL ACCRETION COMPARISON")
    print("=" * 166)
    print(
        f"{'CASE':24s} {'G mdot':>10s} {'TNG/G':>8s} {'BOOST/G':>8s} {'AM/G':>8s} "
        f"{'FFRv/A/G':>9s} {'FFRs/A/G':>9s} {'NEW/G':>8s} "
        f"{'FMACH':>7s} {'FJ':>7s} {'Cconv':>7s} {'Mbulk':>7s} {'Meff':>7s} {'chiJ':>8s}"
    )
    print("-" * 166)
    for r in rows:
        print(
            f"{r['case']:24s} {r['shell_geom_mdot_msun_yr']:10.3e} "
            f"{r['TNG_over_geom']:8.3f} {r['BOOSTED_over_geom']:8.3f} "
            f"{r['AM_over_geom']:8.3f} {r['FFR_VOLUME_over_Aff_geom']:9.3f} "
            f"{r['FFR_SHELL_over_Aff_geom']:9.3f} {r['CONVJ_over_geom']:8.3f} "
            f"{r['FMACH']:7.3f} {r['FJ']:7.3f} {r['CCONV']:7.3f} "
            f"{r['Mbulk']:7.3f} {r['Meff']:7.3f} {r['chiJ']:8.3f}"
        )

    print()
    print(f"Cases with current FFR_CONVJ diagnostics: {len(rows)} / 42")
    print(f"Saved {outfile.name}")
    print("G = unit-efficiency geometric shell FFR; FFRv/A and FFRs/A remove configured A_ff.")
    print("NEW = convergence/J-corrected shell FFR (mode 5).")
    if args.with_reference:
        print("Reference columns were added to the CSV only; they are secondary diagnostics.")
    if missing:
        print()
        print("Cases without current FFR_CONVJ diagnostics:")
        for case in missing:
            print(f"  {case}")


if __name__ == "__main__":
    main()
