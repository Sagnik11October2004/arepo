#!/usr/bin/env python3
"""Generate nested high-resolution reference ICs for resolved BH accretion tests.

These ICs preserve the same physical flow definitions as generate_all_ics.py,
but refine only the central region.  HR4/HR5/HR6 correspond to finest linear
resolutions of 100 pc / (16 * 2^N) for N=4,5,6.

The reference runs use these ICs with a tiny rapidly draining inner absorber.
The absorber rate is *not* treated as a subgrid accretion prescription; the
measured, converged mass reaching it is the numerical reference accretion rate.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
FLOW = HERE.parent
sys.path.insert(0, str(FLOW))
sys.path.insert(0, str(FLOW / "common"))

import flow_ic_utils as U
import generate_all_ics as G


TRUTH_LEVELS = {"HR4": 4, "HR5": 5, "HR6": 6}

# A nested ladder keeps the particle count modest while reaching 0.0977 pc at
# HR6.  The interface radius halves as the linear resolution doubles.
REFINE_RADII_PC = {
    1: 40.0,
    2: 25.0,
    3: 12.5,
    4: 6.25,
    5: 3.125,
    6: 1.5625,
}

BHL_CASES = [("M0p0", 0.0), ("M0p5", 0.5), ("M1p0", 1.0), ("M2p0", 2.0), ("M5p0", 5.0)]
ROT_CASES = [("f0p00", 0.0), ("f0p25", 0.25), ("f0p50", 0.5), ("f0p75", 0.75), ("f1p00", 1.0)]
TURB_CASES = [("M0p5", 0.5), ("M1p0", 1.0), ("M2p0", 2.0), ("M5p0", 5.0)]


def gas_subset(gas: U.GasState, mask: np.ndarray) -> U.GasState:
    return U.GasState(
        gas.coords_kpc[mask].copy(),
        gas.masses_code[mask].copy(),
        gas.velocities_kms[mask].copy(),
        gas.internal_energy_code[mask].copy(),
    )


def gas_concat(a: U.GasState, b: U.GasState) -> U.GasState:
    return U.GasState(
        np.concatenate((a.coords_kpc, b.coords_kpc), axis=0),
        np.concatenate((a.masses_code, b.masses_code), axis=0),
        np.concatenate((a.velocities_kms, b.velocities_kms), axis=0),
        np.concatenate((a.internal_energy_code, b.internal_energy_code), axis=0),
    )


def conserved_totals(gas: U.GasState) -> tuple[float, np.ndarray]:
    mass = float(np.sum(gas.masses_code))
    momentum = np.sum(gas.masses_code[:, None] * gas.velocities_kms, axis=0)
    return mass, momentum


def relative_error(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1.0e-300)


def vector_relative_error(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(a), np.linalg.norm(b), 1.0e-300))


def split_central_region(
    gas: U.GasState,
    parent_dx_kpc: float,
    refine_radius_pc: float,
    density_fn,
    velocity_fn,
    internal_energy_fn,
) -> tuple[U.GasState, dict]:
    """Split every cell centre inside refine_radius_pc into 8 children.

    U.split_state conserves each selected parent's mass and vector momentum
    exactly to floating-point precision.  Unselected cells are carried through
    unchanged.
    """
    rel_pc = U.relative_pc(gas.coords_kpc)
    r_pc = np.linalg.norm(rel_pc, axis=1)
    mask = r_pc < refine_radius_pc
    n_parent = int(np.count_nonzero(mask))
    if n_parent == 0:
        raise RuntimeError(f"no cells selected within {refine_radius_pc} pc")

    before_m, before_p = conserved_totals(gas)

    selected = gas_subset(gas, mask)
    kept = gas_subset(gas, ~mask)
    children = U.split_state(
        selected,
        parent_dx_kpc,
        density_fn,
        velocity_fn,
        internal_energy_fn,
    )
    out = gas_concat(kept, children)

    after_m, after_p = conserved_totals(out)
    report = {
        "refine_radius_pc": float(refine_radius_pc),
        "parents_split": n_parent,
        "cells_before": int(gas.coords_kpc.shape[0]),
        "cells_after": int(out.coords_kpc.shape[0]),
        "mass_relative_error": relative_error(before_m, after_m),
        "momentum_relative_error": vector_relative_error(before_p, after_p),
    }
    return out, report


def make_l0(density_fn, velocity_fn, internal_energy_fn) -> U.GasState:
    coords = U.regular_l0_positions(16)
    rho = density_fn(coords)
    mass = U.masses_from_density_on_l0(coords, rho, 16)
    vel = velocity_fn(coords)
    u = internal_energy_fn(coords)
    return U.GasState(coords, mass, vel, u)


def build_truth_hierarchy(density_fn, velocity_fn, internal_energy_fn) -> tuple[dict[str, U.GasState], list[dict]]:
    gas = make_l0(density_fn, velocity_fn, internal_energy_fn)
    states: dict[str, U.GasState] = {}
    refinement: list[dict] = []

    for level in range(1, 7):
        parent_dx_kpc = U.BOX_KPC / (16.0 * 2.0 ** (level - 1))
        gas, report = split_central_region(
            gas,
            parent_dx_kpc,
            REFINE_RADII_PC[level],
            density_fn,
            velocity_fn,
            internal_energy_fn,
        )
        report["split_level"] = level
        report["child_dx_pc"] = U.BOX_PC / (16.0 * 2.0 ** level)
        refinement.append(report)

        label = f"HR{level}"
        if label in TRUTH_LEVELS:
            states[label] = gas

    return states, refinement


def common_metadata(
    family: str,
    level: str,
    control_name: str,
    control_value: float,
    gas: U.GasState,
    refinement: list[dict],
) -> dict:
    finest_index = TRUTH_LEVELS[level]
    finest_dx_pc = U.BOX_PC / (16.0 * 2.0 ** finest_index)
    net1, in1, n1 = U.shell_flux_msun_per_yr(gas, 10.0, U.R_ACC_PC)
    net2, in2, n2 = U.shell_flux_msun_per_yr(gas, 22.5, U.R_FEEDBACK_PC)

    meta = {
        "family": family,
        "level": level,
        "control_name": control_name,
        "control_value": float(control_value),
        "reference_truth_ic": True,
        "box_pc": U.BOX_PC,
        "bh_mass_msun": U.BH_MASS_MSUN,
        "bh_position_kpc": U.CENTER_KPC.tolist(),
        "bh_velocity_kms": [0.0, 0.0, 0.0],
        "bh_id": int(U.BH_ID),
        "finest_dx_pc": float(finest_dx_pc),
        "recommended_sink_radius_pc": float(4.0 * finest_dx_pc),
        "refinement_ladder_pc": {str(k): float(v) for k, v in REFINE_RADII_PC.items()},
        "refinement_history": refinement[:finest_index],
        "ngas": int(gas.coords_kpc.shape[0]),
        "resolved_flux_diagnostics": {
            "net_10_12p5_msun_per_yr": net1,
            "inward_10_12p5_msun_per_yr": in1,
            "cells_10_12p5": n1,
            "net_22p5_25_msun_per_yr": net2,
            "inward_22p5_25_msun_per_yr": in2,
            "cells_22p5_25": n2,
        },
        "vphi_proxy_within_racc_kms": U.vphi_proxy_kms(gas),
    }
    meta.update(U.basic_state_diagnostics(gas))
    return meta


def write_case(family: str, tag: str, level: str, gas: U.GasState, meta: dict) -> dict:
    stem = f"truth_{family}_{tag}_{level}"
    ic = HERE / family / "ics" / level / f"{stem}.hdf5"
    js = HERE / family / "ics" / level / f"{stem}.metadata.json"
    U.write_arepo_ic(ic, gas, meta)
    U.write_json(js, meta)
    return {
        "family": family,
        "tag": tag,
        "level": level,
        "control_name": meta["control_name"],
        "control_value": meta["control_value"],
        "ic": str(ic.relative_to(FLOW)),
        "metadata": str(js.relative_to(FLOW)),
        "ngas": meta["ngas"],
        "finest_dx_pc": meta["finest_dx_pc"],
        "recommended_sink_radius_pc": meta["recommended_sink_radius_pc"],
    }


def build_bhl(manifest: list[dict]) -> None:
    rho = G.bhl_density_cgs()
    u0 = U.specific_u_from_sound_speed(U.BHL_CS_KMS)

    def density_fn(x):
        return np.full(x.shape[0], rho, dtype=np.float64)

    def u_fn(x):
        return np.full(x.shape[0], u0, dtype=np.float64)

    M_cgs = U.BH_MASS_MSUN * U.M_SUN_CGS

    for tag, mach in BHL_CASES:
        vinf = mach * U.BHL_CS_KMS

        def velocity_fn(x, vinf=vinf):
            v = np.zeros((x.shape[0], 3), dtype=np.float64)
            v[:, 0] = vinf
            return v

        states, refinement = build_truth_hierarchy(density_fn, velocity_fn, u_fn)
        r_bhl_pc = 2.0 * U.G_PC_KMS2_PER_MSUN * U.BH_MASS_MSUN / (U.BHL_CS_KMS**2 + vinf**2)
        rate_cgs = (
            4.0 * np.pi * U.G_CGS**2 * M_cgs**2 * rho
            / (((U.BHL_CS_KMS**2 + vinf**2) ** 1.5) * U.KM_CGS**3)
        )
        rate = rate_cgs * U.YR_CGS / U.M_SUN_CGS

        for level, gas in states.items():
            meta = common_metadata("bhl", level, "mach_infinity", mach, gas, refinement)
            meta.update({
                "sound_speed_kms": U.BHL_CS_KMS,
                "upstream_velocity_kms": vinf,
                "upstream_density_cgs": rho,
                "analytic_bhl_mdot_msun_per_yr": float(rate),
                "r_bhl_pc": float(r_bhl_pc),
                "r_bhl_over_finest_dx": float(r_bhl_pc / meta["finest_dx_pc"]),
                "initial_condition": "uniform gas flowing in +x past a stationary central BH",
            })
            manifest.append(write_case("bhl", tag, level, gas, meta))


def build_rotating(manifest: list[dict]) -> None:
    u0 = U.specific_u_from_temperature(U.SPHERICAL_T_K)

    def u_fn(x):
        return np.full(x.shape[0], u0, dtype=np.float64)

    cs = math.sqrt(U.GAMMA * (U.GAMMA - 1.0) * u0)

    for tag, frot in ROT_CASES:
        density_fn = lambda x, f=frot: U.rotating_density_field(x, f)
        velocity_fn = lambda x, f=frot: U.rotating_velocity_field(x, f)
        states, refinement = build_truth_hierarchy(density_fn, velocity_fn, u_fn)

        for level, gas in states.items():
            meta = common_metadata("rotating", level, "f_rot", frot, gas, refinement)
            meta.update({
                "temperature_K": U.SPHERICAL_T_K,
                "sound_speed_kms": cs,
                "mdot_true_msun_per_yr": U.MDOT_TRUE_MSUN_PER_YR,
                "r_core_pc": U.R_CORE_PC,
                "flow_full_to_pc": U.FLOW_R_FULL_PC,
                "flow_taper_to_pc": U.FLOW_R_TAPER_PC,
                "initial_condition": "spherical supplied inflow plus coherent rotation",
            })
            manifest.append(write_case("rotating", tag, level, gas, meta))


def build_turbulent(manifest: list[dict]) -> None:
    u0 = U.specific_u_from_temperature(U.SPHERICAL_T_K)

    def u_fn(x):
        return np.full(x.shape[0], u0, dtype=np.float64)

    cs = math.sqrt(U.GAMMA * (U.GAMMA - 1.0) * u0)
    density_fn = U.base_spherical_density_field

    l0_coords = U.regular_l0_positions(16)
    l0_mass = U.masses_from_density_on_l0(l0_coords, density_fn(l0_coords), 16)

    for tag, mach in TURB_CASES:
        turb = U.SolenoidalTurbulence(U.TURB_SEED, U.TURB_KMAX)
        turb.calibrate(l0_coords, l0_mass, mach * cs, U.TURB_NORMALIZATION_RADIUS_PC)

        def velocity_fn(x, turb=turb):
            return U.base_spherical_velocity_field(x) + turb(x)

        states, refinement = build_truth_hierarchy(density_fn, velocity_fn, u_fn)

        for level, gas in states.items():
            meta = common_metadata("turbulent", level, "mach_turb", mach, gas, refinement)
            meta.update({
                "temperature_K": U.SPHERICAL_T_K,
                "sound_speed_kms": cs,
                "mdot_true_msun_per_yr": U.MDOT_TRUE_MSUN_PER_YR,
                "turbulence_seed": U.TURB_SEED,
                "target_turbulent_rms_kms": mach * cs,
                "turbulence_normalization_radius_pc": U.TURB_NORMALIZATION_RADIUS_PC,
                "initial_condition": "spherical supplied inflow plus fixed-seed solenoidal turbulence",
            })
            manifest.append(write_case("turbulent", tag, level, gas, meta))


def main() -> None:
    manifest: list[dict] = []
    build_bhl(manifest)
    build_rotating(manifest)
    build_turbulent(manifest)

    path = HERE / "reference_ic_manifest.json"
    path.write_text(json.dumps({"ics": manifest}, indent=2, sort_keys=True) + "\n")

    print(f"Generated {len(manifest)} nested reference ICs")
    for level in TRUTH_LEVELS:
        counts = sorted({x["ngas"] for x in manifest if x["level"] == level})
        dx = U.BOX_PC / (16.0 * 2.0 ** TRUTH_LEVELS[level])
        sink = 4.0 * dx
        print(f"  {level}: finest dx={dx:.8g} pc, recommended sink={sink:.8g} pc, ngas={counts}")
    print(f"Wrote {path.relative_to(FLOW)}")


if __name__ == "__main__":
    main()
