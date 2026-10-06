#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "common"))

from flow_ic_utils import *

BHL_MACHS = [0.0, 0.5, 1.0, 2.0, 5.0]
ROT_FRACTIONS = [0.0, 0.25, 0.5, 0.75, 1.0]
TURB_MACHS = [0.5, 1.0, 2.0, 5.0]


def tag_float(prefix: str, value: float) -> str:
    return f"{prefix}{value:.2f}".replace(".", "p").replace("-", "m").rstrip("0").rstrip("p")


def common_metadata(family: str, level: str, control_name: str, control_value: float, gas: GasState) -> dict:
    net1, in1, n1 = shell_flux_msun_per_yr(gas, 10.0, R_ACC_PC)
    net2, in2, n2 = shell_flux_msun_per_yr(gas, 22.5, R_FEEDBACK_PC)
    d = {
        "family": family,
        "level": level,
        "control_name": control_name,
        "control_value": float(control_value),
        "box_pc": BOX_PC,
        "nside": LEVELS[level],
        "dx_pc": cell_dx_pc(level),
        "bh_mass_msun": BH_MASS_MSUN,
        "bh_position_kpc": CENTER_KPC.tolist(),
        "bh_velocity_kms": [0.0, 0.0, 0.0],
        "bh_id": int(BH_ID),
        "r_acc_pc": R_ACC_PC,
        "r_feedback_pc": R_FEEDBACK_PC,
        "mu": MU,
        "gamma": GAMMA,
        "units": {
            "length": "kpc",
            "mass": "1e10 Msun",
            "velocity": "km/s",
            "internal_energy": "(km/s)^2",
        },
        "resolved_flux_diagnostics": {
            "net_10_12p5_msun_per_yr": net1,
            "inward_10_12p5_msun_per_yr": in1,
            "cells_10_12p5": n1,
            "net_22p5_25_msun_per_yr": net2,
            "inward_22p5_25_msun_per_yr": in2,
            "cells_22p5_25": n2,
            "warning_bhl": "At t=0 a uniform BHL stream has incoming and outgoing hemispheres; these spherical fluxes are geometric diagnostics, not the analytic BHL accretion rate.",
        },
        "vphi_proxy_within_racc_kms": vphi_proxy_kms(gas),
    }
    d.update(basic_state_diagnostics(gas))
    return d


def write_case(root: Path, family: str, level: str, stem: str, gas: GasState, meta: dict) -> tuple[Path, Path]:
    ic = root / family / "ics" / level / f"{stem}.hdf5"
    js = root / family / "ics" / level / f"{stem}.metadata.json"
    write_arepo_ic(ic, gas, meta)
    write_json(js, meta)
    return ic, js


def bhl_density_cgs() -> float:
    mdot = MDOT_TRUE_MSUN_PER_YR * M_SUN_CGS / YR_CGS
    M = BH_MASS_MSUN * M_SUN_CGS
    cs = BHL_CS_KMS * KM_CGS
    v = BHL_FIDUCIAL_MACH * cs
    return mdot * (cs * cs + v * v) ** 1.5 / (4.0 * np.pi * G_CGS**2 * M**2)


def build_bhl(root: Path, manifest: list) -> None:
    rho = bhl_density_cgs()
    nH = X_H * rho / M_P_CGS
    T = temperature_from_sound_speed(BHL_CS_KMS)
    u0 = specific_u_from_sound_speed(BHL_CS_KMS)

    def density_fn(x):
        return np.full(x.shape[0], rho, dtype=np.float64)

    def u_fn(x):
        return np.full(x.shape[0], u0, dtype=np.float64)

    for mach in BHL_MACHS:
        vinf = mach * BHL_CS_KMS
        def velocity_fn(x, vinf=vinf):
            v = np.zeros((x.shape[0], 3), dtype=np.float64)
            v[:, 0] = vinf
            return v
        hierarchy = build_nested_hierarchy(density_fn, velocity_fn, u_fn)
        r_bhl_pc = 2.0 * G_PC_KMS2_PER_MSUN * BH_MASS_MSUN / (BHL_CS_KMS**2 + vinf**2)
        M_cgs = BH_MASS_MSUN * M_SUN_CGS
        rate_cgs = 4.0 * np.pi * G_CGS**2 * M_cgs**2 * rho / (((BHL_CS_KMS**2 + vinf**2) ** 1.5) * KM_CGS**3)
        rate = rate_cgs * YR_CGS / M_SUN_CGS
        tag = f"M{mach:.1f}".replace(".", "p")
        for level, gas in hierarchy.items():
            meta = common_metadata("bhl", level, "mach_infinity", mach, gas)
            meta.update({
                "sound_speed_kms": BHL_CS_KMS,
                "temperature_K": T,
                "upstream_velocity_kms": vinf,
                "upstream_density_cgs": rho,
                "upstream_nH_cm3": nH,
                "analytic_bhl_mdot_msun_per_yr": float(rate),
                "r_bhl_pc": float(r_bhl_pc),
                "r_bhl_over_dx": float(r_bhl_pc / cell_dx_pc(level)),
                "r_bhl_over_racc": float(r_bhl_pc / R_ACC_PC),
                "initial_condition": "uniform gas flowing in +x past a stationary central BH",
            })
            stem = f"bhl_{tag}_{level}"
            ic, js = write_case(root, "bhl", level, stem, gas, meta)
            manifest.append({"ic": str(ic.relative_to(root)), "metadata": str(js.relative_to(root)), **meta})


def build_rotating(root: Path, manifest: list) -> None:
    u0 = specific_u_from_temperature(SPHERICAL_T_K)
    def u_fn(x):
        return np.full(x.shape[0], u0, dtype=np.float64)

    cs = math.sqrt(GAMMA * (GAMMA - 1.0) * u0)
    for frot in ROT_FRACTIONS:
        density_fn = lambda x, f=frot: rotating_density_field(x, f)
        velocity_fn = lambda x, f=frot: rotating_velocity_field(x, f)
        hierarchy = build_nested_hierarchy(density_fn, velocity_fn, u_fn)
        tag = f"f{frot:.2f}".replace(".", "p")
        for level, gas in hierarchy.items():
            vcirc_racc = math.sqrt(G_PC_KMS2_PER_MSUN * BH_MASS_MSUN / R_ACC_PC)
            vphi_racc = frot * vcirc_racc
            meta = common_metadata("rotating", level, "f_rot", frot, gas)
            meta.update({
                "temperature_K": SPHERICAL_T_K,
                "sound_speed_kms": cs,
                "mdot_true_msun_per_yr": MDOT_TRUE_MSUN_PER_YR,
                "r_core_pc": R_CORE_PC,
                "flow_full_to_pc": FLOW_R_FULL_PC,
                "flow_taper_to_pc": FLOW_R_TAPER_PC,
                "ambient_density_fraction": 0.01,
                "velocity_model": "vr=-sqrt((2-frot^2)GM/r_eff), vphi=frot*sqrt(GM/r_eff), tapered from 35 to 40 pc",
                "density_model": "rho=Mdot/(4*pi*r_eff^2*|vr|) inside flow, smooth transition to 1% ambient",
                "vphi_over_cs_at_racc": float(vphi_racc / cs),
                "r_circ_over_racc": float(frot**2),
            })
            stem = f"rot_{tag}_{level}"
            ic, js = write_case(root, "rotating", level, stem, gas, meta)
            manifest.append({"ic": str(ic.relative_to(root)), "metadata": str(js.relative_to(root)), **meta})


def build_turbulent(root: Path, manifest: list) -> None:
    u0 = specific_u_from_temperature(SPHERICAL_T_K)
    cs = math.sqrt(GAMMA * (GAMMA - 1.0) * u0)
    density_fn = base_spherical_density_field
    base_velocity = base_spherical_velocity_field
    def u_fn(x):
        return np.full(x.shape[0], u0, dtype=np.float64)

    # L0 density/masses define the common discrete constraints for removing bulk
    # momentum and coherent rotation. The resulting corrected Fourier field is then
    # evaluated at all child positions; no resolution-specific random field is made.
    l0_coords = regular_l0_positions(16)
    l0_mass = masses_from_density_on_l0(l0_coords, density_fn(l0_coords), 16)

    baseline_hierarchy = build_nested_hierarchy(density_fn, base_velocity, u_fn)

    for mach in TURB_MACHS:
        turb = SolenoidalTurbulence(TURB_SEED, TURB_KMAX)
        turb.calibrate(l0_coords, l0_mass, mach * cs, TURB_NORMALIZATION_RADIUS_PC)
        def velocity_fn(x, turb=turb):
            return base_velocity(x) + turb(x)
        hierarchy = build_nested_hierarchy(density_fn, velocity_fn, u_fn)
        tag = f"M{mach:.1f}".replace(".", "p")
        for level, gas in hierarchy.items():
            baseline = baseline_hierarchy[level]
            rel = relative_pc(gas.coords_kpc)
            r = np.linalg.norm(rel, axis=1)
            mask_norm = r <= TURB_NORMALIZATION_RADIUS_PC
            mask_acc = r <= R_ACC_PC
            dv = gas.velocities_kms - baseline.velocities_kms
            m = gas.masses_code
            actual_rms = math.sqrt(float(np.sum(m[mask_norm] * np.sum(dv[mask_norm]**2, axis=1)) / np.sum(m[mask_norm])))
            rms_acc = math.sqrt(float(np.sum(m[mask_acc] * np.sum(dv[mask_acc]**2, axis=1)) / np.sum(m[mask_acc])))
            meta = common_metadata("turbulent", level, "mach_turb", mach, gas)
            meta.update({
                "temperature_K": SPHERICAL_T_K,
                "sound_speed_kms": cs,
                "mdot_true_msun_per_yr": MDOT_TRUE_MSUN_PER_YR,
                "r_core_pc": R_CORE_PC,
                "flow_full_to_pc": FLOW_R_FULL_PC,
                "flow_taper_to_pc": FLOW_R_TAPER_PC,
                "ambient_density_fraction": 0.01,
                "turbulence_seed": TURB_SEED,
                "turbulence_spectrum": "solenoidal cosine modes with amplitude proportional to |k|^-2, so P_v(k) proportional to k^-4",
                "turbulence_k_integer_max": TURB_KMAX,
                "turbulence_wavelength_range_pc": [BOX_PC / TURB_KMAX, BOX_PC],
                "turbulence_normalization_radius_pc": TURB_NORMALIZATION_RADIUS_PC,
                "target_turbulent_rms_kms": mach * cs,
                "actual_turbulent_rms_within_norm_radius_kms": actual_rms,
                "actual_turbulent_mach_within_norm_radius": actual_rms / cs,
                "actual_turbulent_rms_within_racc_kms": rms_acc,
                "actual_turbulent_mach_within_racc": rms_acc / cs,
                "continuous_field_bulk_subtraction": turb.bulk.tolist(),
                "continuous_field_rotation_omega_per_pc": turb.omega.tolist(),
                "continuous_field_scale_kms": float(turb.scale),
                "mode_table": turb.mode_table(),
            })
            stem = f"turb_{tag}_{level}"
            ic, js = write_case(root, "turbulent", level, stem, gas, meta)
            manifest.append({"ic": str(ic.relative_to(root)), "metadata": str(js.relative_to(root)), **meta})


def main() -> None:
    p = argparse.ArgumentParser(description="Generate the 42 nested AREPO BH accretion flow ICs.")
    p.add_argument("--root", type=Path, default=HERE, help="bh_accretion_flow_tests directory")
    args = p.parse_args()
    root = args.root.resolve()
    manifest = []
    build_bhl(root, manifest)
    build_rotating(root, manifest)
    build_turbulent(root, manifest)
    write_json(root / "ic_manifest.json", {
        "schema_version": 1,
        "count": len(manifest),
        "expected_count": 42,
        "cases": manifest,
    })
    print(f"Generated {len(manifest)} ICs under {root}")
    if len(manifest) != 42:
        raise SystemExit("ERROR: expected exactly 42 ICs")


if __name__ == "__main__":
    main()
