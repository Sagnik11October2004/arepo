#!/usr/bin/env python3
"""Compare the force-corrected shell FFR diagnostic with explicit case references.

Mode 5 itself is evaluated inside AREPO because it uses the local reconstructed
pressure-gradient vector SphP[].Grad.dpress.  This script never derives the
reference from mode 5.

Reference policy
----------------
BHL, gamma=5/3, Mach<1:
    Bondi point-accretor rate with lambda=1/4 from asymptotic IC rho and cs.
BHL, Mach>=1:
    Foglizzo-Ruffert interpolation from asymptotic IC rho, cs, and vinf.
Rotating/turbulent:
    analytic continuum supply used to construct the IC = 1e-4 Msun/yr.
The same-resolution discrete shell flux remains a secondary IC-discretization
diagnostic only.

The script reads the first and last BH_BENCHMARK_ALL lines in each rerun log.
The last sample is the primary comparison reported in the table.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent

G_CGS = 6.67430e-8
M_SUN_CGS = 1.98847e33
YR_CGS = 365.25 * 86400.0
KM_CGS = 1.0e5
LAMBDA_BONDI_GAMMA_5_3 = 0.25
RATE_CODE_TO_MSUN_YR = 10.2271202634


def extract(line: str, key: str) -> float | None:
    m = re.search(
        rf'(?<![A-Za-z0-9_]){re.escape(key)}='
        r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)',
        line,
    )
    return float(m.group(1)) if m else None


def all_lines(log: Path) -> list[str]:
    if not log.exists():
        return []
    out = []
    with log.open("r", errors="replace") as f:
        for line in f:
            if "BH_BENCHMARK_ALL:" in line and "FFR_FORCE=" in line:
                out.append(line)
    return out


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


def stem(family: str, tag: str, level: str) -> str:
    if family == "bhl":
        return f"bhl_{tag}_{level}"
    if family == "rotating":
        return f"rot_{tag}_{level}"
    if family == "turbulent":
        return f"turb_{tag}_{level}"
    raise ValueError(family)


def mdot_tng_physical(mass_msun: float, rho_cgs: float, cs_kms: float) -> float:
    m = mass_msun * M_SUN_CGS
    cs = cs_kms * KM_CGS
    rate_cgs = 4.0 * math.pi * G_CGS**2 * m**2 * rho_cgs / cs**3
    return rate_cgs * YR_CGS / M_SUN_CGS


def bhl_reference(meta: dict) -> dict[str, float | str]:
    mass = float(meta["bh_mass_msun"])
    rho = float(meta["upstream_density_cgs"])
    cs = float(meta["sound_speed_kms"])
    vinf = float(meta["upstream_velocity_kms"])
    mach = vinf / cs
    tng0 = mdot_tng_physical(mass, rho, cs)
    bondi = LAMBDA_BONDI_GAMMA_5_3 * tng0
    if mach < 1.0 - 1.0e-12:
        ref = bondi
        kind = "Bondi-gamma5over3-subsonic"
        status = "strong analytic/numerical point-accretor reference"
    else:
        q = 1.0 + mach * mach
        ref = tng0 * math.sqrt(LAMBDA_BONDI_GAMMA_5_3**2 + mach * mach) / (q * q)
        kind = "Foglizzo-Ruffert-interpolation"
        status = "literature point-accretor baseline; not exact hydrodynamic truth"
    return {
        "reference_mdot_msun_yr": ref,
        "reference_kind": kind,
        "reference_status": status,
        "reference_mach": mach,
        "bondi_gamma5over3_mdot_msun_yr": bondi,
    }


def supply_reference(meta: dict) -> dict[str, float | str]:
    flux = meta["resolved_flux_diagnostics"]
    continuum = float(meta["mdot_true_msun_per_yr"])
    net = float(flux["net_10_12p5_msun_per_yr"])
    inward = float(flux["inward_10_12p5_msun_per_yr"])
    return {
        "reference_mdot_msun_yr": continuum,
        "reference_kind": "analytic-continuum-IC-supply",
        "reference_status": "physical outer supply; not event-horizon BH rate",
        "same_resolution_net_shell_mdot_msun_yr": net,
        "same_resolution_inward_shell_mdot_msun_yr": inward,
        "resolved_net_over_continuum": net / continuum,
    }


def sample(line: str) -> dict[str, float]:
    force = extract(line, "FFR_FORCE")
    geom = extract(line, "FFR_SHELL_GEOM")
    factor = extract(line, "forcefactor")
    if force is None or geom is None:
        raise ValueError("missing FFR_FORCE/FFR_SHELL_GEOM")
    if factor is None:
        factor = force / geom if geom > 0 else 1.0
    return {
        "force_code": force,
        "geom_code": geom,
        "forcefactor": factor,
        "force_mdot_msun_yr": force * RATE_CODE_TO_MSUN_YR,
        "geom_mdot_msun_yr": geom * RATE_CODE_TO_MSUN_YR,
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fields.append(key)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    rows: list[dict] = []
    missing: list[str] = []

    for family, cases in family_cases().items():
        for tag, control in cases:
            for level in ("L0", "L1", "L2"):
                case = stem(family, tag, level)
                meta_path = HERE / family / "ics" / level / f"{case}.metadata.json"
                log = HERE / family / "logs" / "diag" / f"{case}.log"
                meta = json.loads(meta_path.read_text())
                ref = bhl_reference(meta) if family == "bhl" else supply_reference(meta)
                lines = all_lines(log)
                if not lines:
                    missing.append(case)
                    continue

                first = sample(lines[0])
                last = sample(lines[-1])
                mdot_ref = float(ref["reference_mdot_msun_yr"])

                rows.append({
                    "case": case,
                    "family": family,
                    "tag": tag,
                    "control_value": control,
                    "level": level,
                    "samples": len(lines),
                    **ref,
                    "first_force_mdot_msun_yr": first["force_mdot_msun_yr"],
                    "first_force_over_reference": first["force_mdot_msun_yr"] / mdot_ref,
                    "last_force_mdot_msun_yr": last["force_mdot_msun_yr"],
                    "last_force_over_reference": last["force_mdot_msun_yr"] / mdot_ref,
                    "last_geom_mdot_msun_yr": last["geom_mdot_msun_yr"],
                    "last_geom_over_reference": last["geom_mdot_msun_yr"] / mdot_ref,
                    "last_forcefactor": last["forcefactor"],
                    "first_to_last_force_ratio":
                        last["force_mdot_msun_yr"] / first["force_mdot_msun_yr"]
                        if first["force_mdot_msun_yr"] > 0 else math.nan,
                })

    out = HERE / "force_ffr_reference_comparison.csv"
    write_csv(out, rows)

    print()
    print("FORCE-CORRECTED FREE-FALL SHELL / EXPLICIT CASE REFERENCE")
    print("=" * 156)
    print(
        f"{'CASE':24s} {'REF KIND':34s} {'Mdot_ref':>11s} "
        f"{'FORCE/ref':>10s} {'GEOM/ref':>9s} {'factor':>8s} "
        f"{'first/last':>11s} {'Nsamp':>6s}"
    )
    print("-" * 156)
    for r in rows:
        first_last = (
            float(r["first_force_mdot_msun_yr"]) / float(r["last_force_mdot_msun_yr"])
            if float(r["last_force_mdot_msun_yr"]) > 0 else math.nan
        )
        print(
            f"{r['case']:24s} {str(r['reference_kind'])[:34]:34s} "
            f"{float(r['reference_mdot_msun_yr']):11.4e} "
            f"{float(r['last_force_over_reference']):10.4f} "
            f"{float(r['last_geom_over_reference']):9.4f} "
            f"{float(r['last_forcefactor']):8.4f} "
            f"{first_last:11.4f} {int(r['samples']):6d}"
        )

    print()
    print(f"Completed cases with new diagnostics: {len(rows)} / 42")
    print(f"Saved {out.name}")
    if missing:
        print("Missing FFR_FORCE logs:")
        for case in missing:
            print(f"  {case}")
        print("Rebuild, then rerun with FORCE_RERUN=1 ./run_flow_stage.sh diagnostic")

    print()
    print("Reference policy:")
    print("  BHL M<1: gamma=5/3 Bondi rate from asymptotic IC rho and cs.")
    print("  BHL M>=1: Foglizzo-Ruffert interpolation from asymptotic IC rho, cs, vinf.")
    print("  Rotating/turbulent: analytic continuum IC supply = 1e-4 Msun/yr.")
    print("  Same-resolution shell flux is secondary only.")
    print("  No event-horizon truth is assigned to rotating/turbulent cases.")


if __name__ == "__main__":
    main()
