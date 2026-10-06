#!/usr/bin/env python3
"""Evaluate benchmark mode 5 directly from the ICs and compare to explicit references.

Mode 5 is a cell-wise shell/free-fall estimator.  For every gas cell in
0.6 R_acc <= r < R_acc,

    t_ff,i = sqrt(r_i^3 / (G M_central))

    f_i = [1 + (c_s,i^2 + v_perp,i^2 + v_r,out,i^2)
                 / (v_g,i^2 + v_r,in,i^2)]^(-3/2)

    Mdot_eff = C_shell * sum_i m_i f_i / t_ff,i

with
    v_g,i^2   = G M_central / r_i
    v_r,in    = max(-v_r, 0)
    v_r,out   = max(+v_r, 0)
    v_perp^2  = |v_rel|^2 - v_r^2
    C_shell   = sqrt(2) / ln(1/0.6).

The reference is calculated independently of mode 5:
  * BHL, gamma=5/3, Mach<1:
        exact/subsonic Bondi rate with lambda=1/4 from the asymptotic IC state.
  * BHL, Mach>=1:
        Foglizzo-Ruffert interpolation from the asymptotic IC state.
  * Rotating/turbulent:
        the analytic continuum supply used to construct the IC,
        1e-4 Msun/yr.  The same-resolution discrete 10--12.5 pc shell flux
        is reported separately as a discretization diagnostic.

If a newly rerun diagnostic log contains FFR_EFFECTIVE=..., the script also
checks the C implementation against the independent IC calculation.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import re

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent

G_CGS = 6.67430e-8
M_SUN_CGS = 1.98847e33
PC_CGS = 3.0856775814913673e18
YR_CGS = 365.25 * 86400.0
KM_CGS = 1.0e5

G_PC_KMS2_PER_MSUN = G_CGS * M_SUN_CGS / PC_CGS / KM_CGS**2
UNIT_MASS_MSUN = 1.0e10
GAMMA = 5.0 / 3.0
SHELL_INNER_FRACTION = 0.6
SHELL_NORM = math.sqrt(2.0) / math.log(1.0 / SHELL_INNER_FRACTION)
LAMBDA_BONDI_GAMMA_5_3 = 0.25
RATE_CODE_TO_MSUN_YR = 10.2271202634


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


def paths(family: str, tag: str, level: str) -> tuple[Path, Path, Path]:
    s = stem(family, tag, level)
    ic = HERE / family / "ics" / level / f"{s}.hdf5"
    meta = HERE / family / "ics" / level / f"{s}.metadata.json"
    log = HERE / family / "logs" / "diag" / f"{s}.log"
    return ic, meta, log


def periodic_delta(x: np.ndarray, x0: np.ndarray, box: float) -> np.ndarray:
    dx = x - x0[None, :]
    return dx - box * np.rint(dx / box)


def evaluate_ic(ic: Path, meta: dict) -> dict[str, float]:
    if not ic.exists():
        raise FileNotFoundError(ic)

    with h5py.File(ic, "r") as f:
        box_kpc = float(f["Header"].attrs["BoxSize"])
        g = f["PartType0"]
        coords = np.asarray(g["Coordinates"], dtype=np.float64)
        vel = np.asarray(g["Velocities"], dtype=np.float64)
        mass_msun = np.asarray(g["Masses"], dtype=np.float64) * UNIT_MASS_MSUN
        u = np.asarray(g["InternalEnergy"], dtype=np.float64)

        bh = f["PartType5"]
        bh_pos = np.asarray(bh["Coordinates"][0], dtype=np.float64)
        bh_vel = np.asarray(bh["Velocities"][0], dtype=np.float64)
        mcentral_msun = float(bh["Masses"][0]) * UNIT_MASS_MSUN

    rel_pc = periodic_delta(coords, bh_pos, box_kpc) * 1000.0
    r_pc = np.linalg.norm(rel_pc, axis=1)
    racc_pc = float(meta["r_acc_pc"])
    mask = (r_pc >= SHELL_INNER_FRACTION * racc_pc) & (r_pc < racc_pc)
    if not np.any(mask):
        raise RuntimeError(f"{ic}: no shell cells")

    r = r_pc[mask]
    dr = rel_pc[mask]
    dv = vel[mask] - bh_vel[None, :]
    m = mass_msun[mask]
    u_shell = u[mask]

    rhat = dr / r[:, None]
    vr = np.sum(dv * rhat, axis=1)
    dv2 = np.sum(dv * dv, axis=1)
    vperp2 = np.maximum(dv2 - vr * vr, 0.0)

    cs2 = GAMMA * (GAMMA - 1.0) * u_shell
    if np.any(~np.isfinite(cs2)) or np.any(cs2 < 0):
        raise RuntimeError(f"{ic}: invalid sound speed")

    vgrav2 = G_PC_KMS2_PER_MSUN * mcentral_msun / r
    vr_in = np.maximum(-vr, 0.0)
    vr_out = np.maximum(vr, 0.0)
    support2 = cs2 + vperp2 + vr_out**2
    drive2 = vgrav2 + vr_in**2
    factor = (1.0 + support2 / drive2) ** -1.5

    # sqrt(pc^3/[pc (km/s)^2 Msun^-1 * Msun]) = pc/(km/s)
    tff_yr = np.sqrt(r**3 / (G_PC_KMS2_PER_MSUN * mcentral_msun)) * PC_CGS / KM_CGS / YR_CGS
    geom_terms = m / tff_yr
    geom_rate = SHELL_NORM * float(np.sum(geom_terms))
    effective_rate = SHELL_NORM * float(np.sum(geom_terms * factor))
    weighted_factor = effective_rate / geom_rate if geom_rate > 0 else 1.0

    return {
        "shell_cells": int(np.count_nonzero(mask)),
        "central_mass_msun": mcentral_msun,
        "shell_geom_mdot_msun_yr": geom_rate,
        "effective_mdot_msun_yr": effective_rate,
        "effective_support_mean": weighted_factor,
        "factor_min": float(np.min(factor)),
        "factor_median": float(np.median(factor)),
        "factor_max": float(np.max(factor)),
        "mean_vr_kms": float(np.average(vr, weights=m)),
        "mean_abs_vr_kms": float(np.average(np.abs(vr), weights=m)),
        "mean_vperp_kms": float(np.sqrt(np.average(vperp2, weights=m))),
        "mean_cs_kms": float(np.sqrt(np.average(cs2, weights=m))),
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
        "reference_rho_cgs": rho,
        "reference_cs_kms": cs,
        "reference_vinf_kms": vinf,
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
        "reference_status": "physical supply imposed by IC construction; not event-horizon rate",
        "same_resolution_net_shell_mdot_msun_yr": net,
        "same_resolution_inward_shell_mdot_msun_yr": inward,
        "resolved_net_over_continuum": net / continuum,
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

    for family, cases in family_cases().items():
        for tag, control in cases:
            for level in ("L0", "L1", "L2"):
                ic, meta_path, log = paths(family, tag, level)
                meta = json.loads(meta_path.read_text())
                calc = evaluate_ic(ic, meta)
                ref = bhl_reference(meta) if family == "bhl" else supply_reference(meta)
                mdot_ref = float(ref["reference_mdot_msun_yr"])

                row: dict = {
                    "case": stem(family, tag, level),
                    "family": family,
                    "tag": tag,
                    "control_value": control,
                    "level": level,
                    **ref,
                    **calc,
                    "effective_over_reference": calc["effective_mdot_msun_yr"] / mdot_ref,
                    "shell_geom_over_reference": calc["shell_geom_mdot_msun_yr"] / mdot_ref,
                }

                line = first_all_line(log)
                if line is not None:
                    c_eff = extract(line, "FFR_EFFECTIVE")
                    if c_eff is not None:
                        c_eff_msun_yr = c_eff * RATE_CODE_TO_MSUN_YR
                        row["c_effective_mdot_msun_yr"] = c_eff_msun_yr
                        row["c_effective_over_reference"] = c_eff_msun_yr / mdot_ref
                        row["c_over_ic_effective"] = c_eff_msun_yr / calc["effective_mdot_msun_yr"]

                rows.append(row)

    out = HERE / "effective_infall_reference_comparison.csv"
    write_csv(out, rows)

    print()
    print("CELL-WISE EFFECTIVE-INFALL FFR / EXPLICIT CASE REFERENCE")
    print("=" * 158)
    print(
        f"{'CASE':24s} {'REF KIND':34s} {'Mdot_ref':>11s} "
        f"{'EFF/ref':>9s} {'GEOM/ref':>9s} {'<f_eff>':>9s} "
        f"{'f50':>7s} {'vr':>8s} {'vperp':>8s}"
    )
    print("-" * 158)
    for r in rows:
        print(
            f"{r['case']:24s} {str(r['reference_kind'])[:34]:34s} "
            f"{float(r['reference_mdot_msun_yr']):11.4e} "
            f"{float(r['effective_over_reference']):9.4f} "
            f"{float(r['shell_geom_over_reference']):9.4f} "
            f"{float(r['effective_support_mean']):9.4f} "
            f"{float(r['factor_median']):7.4f} "
            f"{float(r['mean_vr_kms']):8.3f} "
            f"{float(r['mean_vperp_kms']):8.3f}"
        )

    print()
    c_rows = [r for r in rows if "c_over_ic_effective" in r]
    print(f"Completed cases: {len(rows)}")
    print(f"Saved {out.name}")
    if c_rows:
        worst = max(abs(float(r["c_over_ic_effective"]) - 1.0) for r in c_rows)
        print(f"C/log cross-check cases: {len(c_rows)}, worst |C/IC-1| = {worst:.3e}")
    else:
        print("C/log cross-check: no FFR_EFFECTIVE diagnostics yet; rerun diagnostics after rebuilding.")

    print()
    print("Reference policy:")
    print("  BHL M<1: gamma=5/3 Bondi rate from asymptotic IC rho and cs.")
    print("  BHL M>=1: Foglizzo-Ruffert interpolation from asymptotic IC rho, cs, vinf.")
    print("  Rotating/turbulent: analytic continuum IC supply = 1e-4 Msun/yr.")
    print("  Same-resolution shell flux is secondary only; it measures IC discretization.")
    print("  No event-horizon truth is assigned to rotating/turbulent cases.")


if __name__ == "__main__":
    main()
