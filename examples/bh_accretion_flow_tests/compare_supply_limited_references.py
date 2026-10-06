#!/usr/bin/env python3
"""Reconstruct the supply-limited hybrid from completed diagnostic logs.

Mode 5:
    HYBRID_FFR = min(FFR_SHELL_GEOM, CAPTURE_CEILING)

For this gamma=5/3 benchmark suite:
  * Mach < 1: CAPTURE_CEILING = exact/subsonic Bondi rate = 0.25*TNG.
  * Mach >= 1: Ruffert-Arnett interpolation
        TNG * sqrt(lambda^2 + M^2) / (1 + M^2)^2, lambda=0.25.

For rotating/turbulent cases, the reference column is the resolved L2 net
supply through the 10--12.5 pc shell. This is a capture/supply reference, not
a claim about the event-horizon accretion rate.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
RATE_CODE_TO_MSUN_YR = 10.2271202634
A_FF_OLD_LOGS = 1.0e-3
LAMBDA_BONDI_GAMMA_5_3 = 0.25

BASE_KEYS = ("TNG", "BOOSTED", "AM", "FFR_VOLUME", "FFR_SHELL")


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


def metadata(family: str, tag: str, level: str) -> dict:
    path = HERE / family / "ics" / level / f"{stem(family, tag, level)}.metadata.json"
    return json.loads(path.read_text())


def family_cases() -> dict[str, list[tuple[str, float]]]:
    return {
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


def l2_supply_reference(family: str, tag: str) -> float:
    m = metadata(family, tag, "L2")
    return float(m["resolved_flux_diagnostics"]["net_10_12p5_msun_per_yr"])


def capture_ceiling_code(tng_code: float, cs: float, vrel: float) -> tuple[float, float, str]:
    if not (cs > 0) or vrel < 0:
        raise ValueError(f"invalid cs/vrel: {cs}, {vrel}")
    mach = vrel / cs
    if mach < 1.0:
        return (
            LAMBDA_BONDI_GAMMA_5_3 * tng_code,
            mach,
            "exact/subsonic-Bondi-gamma5over3",
        )

    q = 1.0 + mach * mach
    rate = (
        tng_code
        * math.sqrt(LAMBDA_BONDI_GAMMA_5_3 ** 2 + mach * mach)
        / (q * q)
    )
    return rate, mach, "Ruffert-Arnett-gamma5over3"


def rates_from_log(line: str) -> dict[str, float | str]:
    code: dict[str, float] = {}
    for key in BASE_KEYS:
        x = extract(line, key)
        if x is None:
            raise ValueError(f"missing {key}=... in BH_BENCHMARK_ALL")
        code[key] = x

    cs = extract(line, "cs")
    vrel = extract(line, "vrel")
    if cs is None or vrel is None:
        raise ValueError("missing cs/vrel")

    shell_geom = extract(line, "FFR_SHELL_GEOM")
    if shell_geom is None:
        # The completed pre-hybrid diagnostic matrix used A_ff=1e-3.
        shell_geom = code["FFR_SHELL"] / A_FF_OLD_LOGS

    capture_code, mach, ceiling_kind = capture_ceiling_code(code["TNG"], cs, vrel)
    hybrid_code = min(shell_geom, capture_code)
    supply_limiter = hybrid_code / shell_geom if shell_geom > 0 else 1.0
    branch = "supply" if shell_geom <= capture_code else "capture"

    out: dict[str, float | str] = {
        **code,
        "FFR_SHELL_GEOM": shell_geom,
        "CAPTURE_CEILING": capture_code,
        "HYBRID_FFR": hybrid_code,
        "bulk_mach": mach,
        "supply_limiter": supply_limiter,
        "limiting_branch": branch,
        "ceiling_kind": ceiling_kind,
    }
    return out


def reference_for_case(family: str, tag: str, rates: dict[str, float | str]) -> dict[str, float | str]:
    if family == "bhl":
        ref = float(rates["CAPTURE_CEILING"]) * RATE_CODE_TO_MSUN_YR
        return {
            "reference_mdot_msun_yr": ref,
            "reference_kind": str(rates["ceiling_kind"]),
            "reference_status": (
                "strong theory reference"
                if float(rates["bulk_mach"]) < 1.0
                else "literature interpolation baseline; not exact hydrodynamic truth"
            ),
        }

    supply = l2_supply_reference(family, tag)
    return {
        "reference_mdot_msun_yr": supply,
        "reference_kind": "resolved-L2-supply-at-Racc",
        "reference_status": "capture/supply reference, not BH-horizon truth",
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields: list[str] = []
    seen: set[str] = set()
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
    rows: list[dict] = []

    for family, cases in family_cases().items():
        for tag, value in cases:
            for level in ("L0", "L1", "L2"):
                case = stem(family, tag, level)
                log = HERE / family / "logs" / "diag" / f"{case}.log"
                line = first_all_line(log)
                if line is None:
                    continue

                rates = rates_from_log(line)
                ref = reference_for_case(family, tag, rates)
                mdot_ref = float(ref["reference_mdot_msun_yr"])

                row: dict = {
                    "case": case,
                    "family": family,
                    "tag": tag,
                    "control_value": value,
                    "level": level,
                    "bulk_mach": rates["bulk_mach"],
                    "limiting_branch": rates["limiting_branch"],
                    "supply_limiter": rates["supply_limiter"],
                    **ref,
                }

                for key in (*BASE_KEYS, "FFR_SHELL_GEOM", "CAPTURE_CEILING", "HYBRID_FFR"):
                    mdot = float(rates[key]) * RATE_CODE_TO_MSUN_YR
                    row[f"{key}_mdot_msun_yr"] = mdot
                    row[f"{key}_over_reference"] = mdot / mdot_ref

                rows.append(row)

    write_csv(HERE / "supply_limited_reference_comparison.csv", rows)

    print()
    print("SUPPLY-LIMITED HYBRID / REFERENCE")
    print("=" * 132)
    print(
        f"{'CASE':24s} {'REF KIND':31s} {'Mdot_ref':>11s} "
        f"{'HYBRID/ref':>11s} {'branch':>9s} {'Mbulk':>8s} {'supplylim':>10s}"
    )
    print("-" * 132)
    for r in rows:
        print(
            f"{r['case']:24s} {str(r['reference_kind'])[:31]:31s} "
            f"{float(r['reference_mdot_msun_yr']):11.4e} "
            f"{float(r['HYBRID_FFR_over_reference']):11.4f} "
            f"{str(r['limiting_branch']):>9s} "
            f"{float(r['bulk_mach']):8.4f} "
            f"{float(r['supply_limiter']):10.4f}"
        )

    print()
    print(f"Completed cases: {len(rows)}")
    print("Saved supply_limited_reference_comparison.csv")
    print()
    print("Interpretation:")
    print("  BHL: hybrid is capped by the gamma=5/3 point-accretor ceiling.")
    print("  Rotating/turbulent: comparison reference is resolved supply at Racc.")
    print("  Angular-momentum transport is NOT included in shell FFR or this hybrid.")


if __name__ == "__main__":
    main()
