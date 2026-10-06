#!/usr/bin/env python3
"""Compare the six raw accretion estimators with defensible reference rates.

Reference hierarchy:
  * BHL, Mach < 1, gamma=5/3: exact/subsonic Bondi rate = 0.25*TNG-Bondi.
  * BHL, Mach >= 1: standard BHL interpolation is used as a literature
    baseline, not labelled exact truth.
  * rotating/turbulent: the L2 resolved net shell flux at 10--12.5 pc is used
    as the reference *supply* rate.  It is not claimed to be the horizon rate.

MACH_FFR is read from new logs when available.  For old completed logs it is
reconstructed exactly from FFR_SHELL, BHFreeFallA=1e-3, cs and vrel.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
RATE_CODE_TO_MSUN_YR = 10.2271202634
A_FF = 1.0e-3

MODE_KEYS = ("TNG", "BOOSTED", "AM", "FFR_VOLUME", "FFR_SHELL", "MACH_FFR")


def extract(line: str, key: str) -> float | None:
    m = re.search(
        rf'(?<![A-Za-z0-9_]){re.escape(key)}='
        r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)',
        line,
    )
    return float(m.group(1)) if m else None


def first_all_line(log: Path) -> str | None:
    if not log.exists():
        return None
    with log.open("r", errors="replace") as f:
        for line in f:
            if "BH_BENCHMARK_ALL:" in line:
                return line
    return None


def stem(family: str, tag: str, level: str) -> str:
    if family == "bhl":
        return f"bhl_{tag}_{level}"
    if family == "rotating":
        return f"rot_{tag}_{level}"
    if family == "turbulent":
        return f"turb_{tag}_{level}"
    raise ValueError(family)


def metadata_path(family: str, tag: str, level: str) -> Path:
    return HERE / family / "ics" / level / f"{stem(family, tag, level)}.metadata.json"


def metadata(family: str, tag: str, level: str) -> dict:
    return json.loads(metadata_path(family, tag, level).read_text())


def family_cases() -> dict[str, list[tuple[str, float]]]:
    return {
        "bhl": [
            ("M0p0", 0.0), ("M0p5", 0.5), ("M1p0", 1.0),
            ("M2p0", 2.0), ("M5p0", 5.0),
        ],
        "rotating": [
            ("f0p00", 0.0), ("f0p25", 0.25), ("f0p50", 0.5),
            ("f0p75", 0.75), ("f1p00", 1.0),
        ],
        "turbulent": [
            ("M0p5", 0.5), ("M1p0", 1.0), ("M2p0", 2.0), ("M5p0", 5.0),
        ],
    }


def l2_supply_reference(family: str, tag: str) -> float:
    m = metadata(family, tag, "L2")
    return float(m["resolved_flux_diagnostics"]["net_10_12p5_msun_per_yr"])


def rates_from_log(line: str) -> tuple[dict[str, float], float]:
    code = {}
    for key in ("TNG", "BOOSTED", "AM", "FFR_VOLUME", "FFR_SHELL"):
        x = extract(line, key)
        if x is None:
            raise ValueError(f"missing {key}=... in BH_BENCHMARK_ALL")
        code[key] = x

    cs = extract(line, "cs")
    vrel = extract(line, "vrel")
    if cs is None or vrel is None or not (cs > 0) or vrel < 0:
        raise ValueError("missing/invalid cs or vrel")

    machsup = (1.0 + (vrel / cs) ** 2) ** -1.5
    mffr_logged = extract(line, "MACH_FFR")
    code["MACH_FFR"] = (
        mffr_logged
        if mffr_logged is not None
        else (code["FFR_SHELL"] / A_FF) * machsup
    )

    return {k: v * RATE_CODE_TO_MSUN_YR for k, v in code.items()}, machsup


def reference_for_case(family: str, tag: str, value: float, rates: dict[str, float]) -> dict:
    if family == "bhl":
        if value < 1.0:
            # For gamma=5/3, resolved subsonic wind accretion tends to the
            # spherical Bondi rate and is independent of subsonic wind speed.
            ref = 0.25 * rates["TNG"]
            return {
                "reference_mdot_msun_yr": ref,
                "reference_kind": "exact/subsonic-Bondi-gamma5over3",
                "reference_status": "strong-theory-reference",
                "reference_low_msun_yr": ref,
                "reference_high_msun_yr": ref,
            }

        # At trans/supersonic Mach there is no universal exact algebraic rate.
        # The standard BHL interpolation is retained as the central literature
        # baseline.  Do not relabel it exact truth; resolved flows can differ.
        m = metadata("bhl", tag, "L2")
        ref = float(m["analytic_bhl_mdot_msun_per_yr"])
        if value >= 3.0:
            lo, hi = 0.1 * ref, 1.3 * ref
            status = "high-Mach-baseline; resolved instabilities can strongly reduce rate"
        else:
            lo, hi = 0.7 * ref, 1.3 * ref
            status = "BHL-literature-baseline; ~30pct modelling uncertainty"
        return {
            "reference_mdot_msun_yr": ref,
            "reference_kind": "standard-BHL-baseline",
            "reference_status": status,
            "reference_low_msun_yr": lo,
            "reference_high_msun_yr": hi,
        }

    supply = l2_supply_reference(family, tag)
    return {
        "reference_mdot_msun_yr": supply,
        "reference_kind": "resolved-L2-supply-at-Racc",
        "reference_status": "capture/supply reference, not BH-horizon truth",
        "reference_low_msun_yr": supply,
        "reference_high_msun_yr": supply,
    }


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
    rows = []

    for family, cases in family_cases().items():
        for tag, value in cases:
            for level in ("L0", "L1", "L2"):
                case = stem(family, tag, level)
                log = HERE / family / "logs" / "diag" / f"{case}.log"
                line = first_all_line(log)
                if line is None:
                    continue

                rates, machsup = rates_from_log(line)
                ref = reference_for_case(family, tag, value, rates)
                mdot_ref = ref["reference_mdot_msun_yr"]

                row = {
                    "case": case,
                    "family": family,
                    "tag": tag,
                    "control_value": value,
                    "level": level,
                    "bulk_mach_suppressor": machsup,
                    **ref,
                }
                for mode in MODE_KEYS:
                    row[f"{mode}_mdot_msun_yr"] = rates[mode]
                    row[f"{mode}_over_reference"] = rates[mode] / mdot_ref

                rows.append(row)

    write_csv(HERE / "mach_ffr_reference_comparison.csv", rows)

    print()
    print("MACH-SUPPRESSED SHELL FFR / REFERENCE")
    print("=" * 116)
    print(f"{'CASE':24s} {'REF KIND':31s} {'Mdot_ref':>11s} {'MACH_FFR/ref':>14s} {'machsup':>10s}")
    print("-" * 116)
    for r in rows:
        print(
            f"{r['case']:24s} {r['reference_kind'][:31]:31s} "
            f"{r['reference_mdot_msun_yr']:11.4e} "
            f"{r['MACH_FFR_over_reference']:14.4f} "
            f"{r['bulk_mach_suppressor']:10.4f}"
        )

    print()
    print(f"Completed cases: {len(rows)}")
    print("Saved mach_ffr_reference_comparison.csv")
    print()
    print("Interpretation:")
    print("  BHL M<1: reference is the gamma=5/3 subsonic Bondi rate.")
    print("  BHL M>=1: reference is a literature BHL baseline, not exact truth.")
    print("  Rotating/turbulent: reference is resolved supply at Racc, not horizon accretion.")


if __name__ == "__main__":
    main()
