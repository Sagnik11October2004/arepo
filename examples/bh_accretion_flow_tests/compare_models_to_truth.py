#!/usr/bin/env python3
"""Compare the five coarse subgrid estimators with converged reference rates."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import re


FLOW = Path(__file__).resolve().parent
TRUTH = FLOW / "reference_truth" / "reference_truth_summary.json"

RATE_CODE_TO_MSUN_YR = 10.2271202634
FFR_A = 1.0e-3
MODES = ("TNG", "BOOSTED", "AM", "FFR_VOLUME", "FFR_SHELL")


def first_all_line(log: Path) -> str | None:
    if not log.exists():
        return None
    with log.open("r", errors="replace") as f:
        for line in f:
            if "BH_BENCHMARK_ALL:" in line:
                return line
    return None


def extract(line: str, key: str) -> float:
    m = re.search(
        rf'(?<![A-Za-z0-9_]){re.escape(key)}='
        r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)',
        line,
    )
    if not m:
        raise ValueError(f"missing {key}=... in benchmark line")
    return float(m.group(1))


def case_stem(family: str, tag: str, level: str) -> str:
    if family == "bhl":
        return f"bhl_{tag}_{level}"
    if family == "rotating":
        return f"rot_{tag}_{level}"
    if family == "turbulent":
        return f"turb_{tag}_{level}"
    raise ValueError(family)


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
    if not TRUTH.exists():
        raise SystemExit(
            "Missing converged-reference summary. Run "
            "python3 reference_truth/analyze_reference_truth.py first."
        )

    summary = json.loads(TRUTH.read_text())["cases"]
    truth = {
        (x["family"], x["tag"]): x
        for x in summary
        if x.get("actual_mdot_msun_yr") is not None
    }

    rows = []

    for (family, tag), ref in sorted(truth.items()):
        actual = float(ref["actual_mdot_msun_yr"])

        for level in ("L0", "L1", "L2"):
            stem = case_stem(family, tag, level)
            log = FLOW / family / "logs" / "diag" / f"{stem}.log"
            line = first_all_line(log)
            if line is None:
                continue

            rates = {
                mode: extract(line, mode) * RATE_CODE_TO_MSUN_YR
                for mode in MODES
            }

            row = {
                "family": family,
                "tag": tag,
                "level": level,
                "actual_mdot_msun_yr": actual,
                "truth_HR5_to_HR6_relative_change":
                    ref.get("HR5_to_HR6_relative_change"),
                "truth_HR6_variability_fraction":
                    ref.get("HR6_variability_fraction"),
            }

            for mode in MODES:
                row[f"{mode}_mdot_msun_yr"] = rates[mode]
                row[f"{mode}_over_actual"] = rates[mode] / actual

            # These two extra columns diagnose FFR geometry independent of its
            # explicit A_ff=1e-3 efficiency.  They are not the configured model
            # rates; the unscaled FFR_*_over_actual columns above are.
            row["FFR_VOLUME_over_Aff_actual"] = rates["FFR_VOLUME"] / (FFR_A * actual)
            row["FFR_SHELL_over_Aff_actual"] = rates["FFR_SHELL"] / (FFR_A * actual)

            rows.append(row)

    rows.sort(key=lambda r: (r["family"], r["tag"], r["level"]))
    write_csv(FLOW / "models_vs_actual.csv", rows)

    print()
    print("SUBGRID ESTIMATORS / CONVERGED HIGH-RESOLUTION ACTUAL RATE")
    print("=" * 126)
    print(
        f"{'CASE':24s} {'actual':>11s} {'TNG':>9s} {'BOOST':>9s} {'AM':>9s} "
        f"{'FFRv':>9s} {'FFRs':>9s}"
    )
    print("-" * 126)

    for r in rows:
        case = f"{case_stem(r['family'], r['tag'], r['level'])}"
        print(
            f"{case:24s}"
            f" {r['actual_mdot_msun_yr']:11.4e}"
            f" {r['TNG_over_actual']:9.3f}"
            f" {r['BOOSTED_over_actual']:9.3f}"
            f" {r['AM_over_actual']:9.3f}"
            f" {r['FFR_VOLUME_over_actual']:9.3f}"
            f" {r['FFR_SHELL_over_actual']:9.3f}"
        )

    print()
    print(f"Rows with both a converged truth rate and completed coarse diagnostic: {len(rows)}")
    print("Saved models_vs_actual.csv")
    print("FFRv/FFRs above are the configured rates including A_ff=1e-3.")
    print("The CSV also contains FFR/(A_ff * actual) for geometry-only diagnosis.")


if __name__ == "__main__":
    main()
