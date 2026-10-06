#!/usr/bin/env python3
"""Compare the unified shell-support mode with case-appropriate references.

Mode 5 reconstructed from existing diagnostic logs:
    M_dyn = M_central + M_gas(<R_acc)
    Mdot_shell,dyn = Mdot_shell,geom * sqrt(M_dyn/M_central)
    v_g^2 = G M_dyn/R_acc
    f_support = [1 + (c_s^2 + v_bulk^2)/v_g^2]^(-3/2)
    Mdot_support = Mdot_shell,dyn * f_support

Reference hierarchy:
  BHL, gamma=5/3, Mach<1:
      spherical Bondi rate with lambda=1/4, computed from the asymptotic IC
      density and sound speed. This is the strongest analytic/numerical
      reference in the suite.
  BHL, Mach>=1:
      Foglizzo-Ruffert/Ruffert-Arnett-style interpolation, evaluated from the
      asymptotic IC state. It is a literature baseline, not exact truth.
  Rotating/turbulent:
      the actual discrete net inward shell flux of THAT resolution through
      10--12.5 pc. The continuum IC target (1e-4 Msun/yr) is also reported
      separately. These are capture/supply references, not event-horizon rates.

No analytic event-horizon "truth" is asserted for rotating/turbulent runs,
because angular-momentum transport and unresolved disc physics determine the
final BH growth rate.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent

# Physical constants consistent with the IC generator.
G_CGS = 6.67430e-8
M_SUN_CGS = 1.98847e33
PC_CGS = 3.0856775814913673e18
YR_CGS = 365.25 * 86400.0
KM_CGS = 1.0e5
G_PC_KMS2_PER_MSUN = G_CGS * M_SUN_CGS / PC_CGS / (KM_CGS**2)

UNIT_MASS_MSUN = 1.0e10
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
    standard_bhl = tng0 / (1.0 + mach * mach) ** 1.5

    if mach < 1.0 - 1.0e-12:
        ref = bondi
        kind = "Bondi-gamma5over3-subsonic"
        status = (
            "strong analytic/numerical point-accretor reference; "
            "subsonic rate is independent of wind Mach for gamma=5/3"
        )
        lo = hi = ref
    else:
        q = 1.0 + mach * mach
        ref = tng0 * math.sqrt(LAMBDA_BONDI_GAMMA_5_3**2 + mach * mach) / (q * q)
        kind = "Foglizzo-Ruffert-interpolation"
        status = (
            "literature point-accretor baseline, not exact hydrodynamic truth; "
            "resolved shocks/instabilities can shift the mean rate"
        )
        # Use a transparent broad comparison band, not a statistical error bar.
        # Around classical resolved tests the interpolation is commonly an
        # order-unity (~30%) guide; high-Mach adiabatic instabilities can lower
        # the rate much more, so widen the lower side for M>=3.
        if mach >= 3.0:
            lo, hi = 0.1 * ref, 1.3 * ref
        else:
            lo, hi = 0.7 * ref, 1.3 * ref

    return {
        "reference_mdot_msun_yr": ref,
        "reference_low_msun_yr": lo,
        "reference_high_msun_yr": hi,
        "reference_kind": kind,
        "reference_status": status,
        "bondi_gamma5over3_mdot_msun_yr": bondi,
        "standard_bhl_mdot_msun_yr": standard_bhl,
    }


def supply_reference(meta: dict) -> dict[str, float | str]:
    flux = meta["resolved_flux_diagnostics"]
    net = float(flux["net_10_12p5_msun_per_yr"])
    inward = float(flux["inward_10_12p5_msun_per_yr"])
    continuum = float(meta["mdot_true_msun_per_yr"])
    return {
        "reference_mdot_msun_yr": net,
        "reference_low_msun_yr": net,
        "reference_high_msun_yr": net,
        "reference_kind": "same-resolution-resolved-shell-supply",
        "reference_status": (
            "direct capture/supply reference at Racc; not an event-horizon BH rate"
        ),
        "continuum_target_mdot_msun_yr": continuum,
        "resolved_inward_shell_mdot_msun_yr": inward,
        "resolved_net_over_continuum": net / continuum,
    }


def reconstruct_mode5(line: str, meta: dict) -> dict[str, float | str]:
    code: dict[str, float] = {}
    for key in BASE_KEYS:
        x = extract(line, key)
        if x is None:
            raise ValueError(f"missing {key}=... in BH_BENCHMARK_ALL")
        code[key] = x

    cs = extract(line, "cs")
    vrel = extract(line, "vrel")
    mgas_code = extract(line, "Mgas")
    if cs is None or vrel is None or mgas_code is None:
        raise ValueError("missing cs/vrel/Mgas")

    shell_geom = extract(line, "FFR_SHELL_GEOM")
    if shell_geom is None:
        # Pre-mode-5 completed logs used A_ff=1e-3.
        shell_geom = code["FFR_SHELL"] / A_FF_OLD_LOGS

    mcentral_msun = float(meta["bh_mass_msun"])
    mgas_msun = mgas_code * UNIT_MASS_MSUN
    mdyn_msun = mcentral_msun + mgas_msun
    racc_pc = float(meta["r_acc_pc"])

    mass_scale = math.sqrt(mdyn_msun / mcentral_msun)
    shell_dyn = shell_geom * mass_scale
    vg2 = G_PC_KMS2_PER_MSUN * mdyn_msun / racc_pc
    vg = math.sqrt(vg2)
    support = (1.0 + (cs * cs + vrel * vrel) / vg2) ** -1.5
    support_ffr = shell_dyn * support

    return {
        **code,
        "FFR_SHELL_GEOM": shell_geom,
        "SHELL_DYN": shell_dyn,
        "SHELL_SUPPORT": support_ffr,
        "support_factor": support,
        "gravity_speed_kms": vg,
        "central_mass_msun": mcentral_msun,
        "gas_mass_within_racc_msun": mgas_msun,
        "dynamical_mass_msun": mdyn_msun,
        "bulk_mach": vrel / cs if cs > 0 else math.inf,
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

                meta = metadata(family, tag, level)
                rates = reconstruct_mode5(line, meta)
                ref = bhl_reference(meta) if family == "bhl" else supply_reference(meta)
                mdot_ref = float(ref["reference_mdot_msun_yr"])

                row: dict = {
                    "case": case,
                    "family": family,
                    "tag": tag,
                    "control_value": value,
                    "level": level,
                    "bulk_mach": rates["bulk_mach"],
                    "central_mass_msun": rates["central_mass_msun"],
                    "gas_mass_within_racc_msun": rates["gas_mass_within_racc_msun"],
                    "dynamical_mass_msun": rates["dynamical_mass_msun"],
                    "gravity_speed_kms": rates["gravity_speed_kms"],
                    "support_factor": rates["support_factor"],
                    **ref,
                }

                for key in (*BASE_KEYS, "FFR_SHELL_GEOM", "SHELL_DYN", "SHELL_SUPPORT"):
                    mdot = float(rates[key]) * RATE_CODE_TO_MSUN_YR
                    row[f"{key}_mdot_msun_yr"] = mdot
                    row[f"{key}_over_reference"] = mdot / mdot_ref if mdot_ref != 0 else math.nan

                rows.append(row)

    out = HERE / "shell_support_reference_comparison.csv"
    write_csv(out, rows)

    print()
    print("UNIFIED SHELL-SUPPORT FFR / CASE-APPROPRIATE REFERENCE")
    print("=" * 150)
    print(
        f"{'CASE':24s} {'REF KIND':34s} {'Mdot_ref':>11s} "
        f"{'SUPPORT/ref':>12s} {'support':>9s} {'vg[km/s]':>9s} "
        f"{'Mdyn[Msun]':>12s}"
    )
    print("-" * 150)
    for r in rows:
        print(
            f"{r['case']:24s} {str(r['reference_kind'])[:34]:34s} "
            f"{float(r['reference_mdot_msun_yr']):11.4e} "
            f"{float(r['SHELL_SUPPORT_over_reference']):12.4f} "
            f"{float(r['support_factor']):9.4f} "
            f"{float(r['gravity_speed_kms']):9.3f} "
            f"{float(r['dynamical_mass_msun']):12.4e}"
        )

    print()
    print(f"Completed cases: {len(rows)}")
    print(f"Saved {out.name}")
    print()
    print("Reference policy:")
    print("  BHL M<1: gamma=5/3 Bondi rate from asymptotic IC state.")
    print("  BHL M>=1: Foglizzo-Ruffert interpolation baseline; not exact truth.")
    print("  Rotating/turbulent: same-resolution resolved net shell supply at Racc.")
    print("  Continuum 1e-4 Msun/yr target is retained separately for rotating/turbulent.")
    print("  No event-horizon truth is assigned to rotating/turbulent cases.")


if __name__ == "__main__":
    main()
