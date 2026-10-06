#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "common"))
from flow_ic_utils import *

EXPECTED = {
    "bhl": [(f"M{x:.1f}".replace(".", "p"), x) for x in (0.0, 0.5, 1.0, 2.0, 5.0)],
    "rotating": [(f"f{x:.2f}".replace(".", "p"), x) for x in (0.0, 0.25, 0.5, 0.75, 1.0)],
    "turbulent": [(f"M{x:.1f}".replace(".", "p"), x) for x in (0.5, 1.0, 2.0, 5.0)],
}


def load_ic(path: Path) -> tuple[GasState, dict]:
    with h5py.File(path, "r") as f:
        h = f["Header"].attrs
        n = np.array(h["NumPart_Total"], dtype=np.uint64)
        if int(n[5]) != 1 or np.sum(n[1:5]) != 0:
            raise AssertionError(f"{path}: expected exactly one Type-5 BH and no other collisionless types")
        for dset in ("PartType0/Coordinates", "PartType0/Masses", "PartType0/Velocities", "PartType0/InternalEnergy",
                     "PartType5/Coordinates", "PartType5/Masses", "PartType5/Velocities"):
            if f[dset].dtype != np.dtype("float64"):
                raise AssertionError(f"{path}: {dset} is {f[dset].dtype}, expected float64")
        gas_ids = np.array(f["PartType0/ParticleIDs"], dtype=np.uint64)
        if gas_ids.shape != (int(n[0]),) or not np.array_equal(gas_ids, np.arange(1, int(n[0])+1, dtype=np.uint64)):
            raise AssertionError(f"{path}: gas IDs are not deterministic sequential IDs")
        gas = GasState(
            np.array(f["PartType0/Coordinates"], dtype=np.float64),
            np.array(f["PartType0/Masses"], dtype=np.float64),
            np.array(f["PartType0/Velocities"], dtype=np.float64),
            np.array(f["PartType0/InternalEnergy"], dtype=np.float64),
        )
        if gas.coords_kpc.shape[0] != int(n[0]):
            raise AssertionError(f"{path}: header gas count mismatch")
        bh_pos = np.array(f["PartType5/Coordinates"], dtype=np.float64)
        bh_vel = np.array(f["PartType5/Velocities"], dtype=np.float64)
        bh_id = np.array(f["PartType5/ParticleIDs"], dtype=np.uint64)
        bh_mass = np.array(f["PartType5/Masses"], dtype=np.float64)
    with path.with_suffix(".metadata.json").open("r", encoding="utf-8") as q:
        meta = json.load(q)
    if meta.get("level") not in LEVELS:
        raise AssertionError(f"{path}: invalid/missing metadata level")
    if not np.array_equal(bh_id, np.array([BH_ID], dtype=np.uint64)):
        raise AssertionError(f"{path}: wrong BH ID")
    if not np.allclose(bh_pos, CENTER_KPC.reshape(1,3), rtol=0, atol=1e-15):
        raise AssertionError(f"{path}: wrong BH position")
    if not np.allclose(bh_vel, 0.0, rtol=0, atol=0):
        raise AssertionError(f"{path}: nonzero BH velocity")
    if not np.allclose(bh_mass, BH_MASS_CODE, rtol=1e-15, atol=0):
        raise AssertionError(f"{path}: wrong BH mass")
    return gas, meta


def finite_positive_checks(path: Path, gas: GasState, level: str) -> None:
    if gas.coords_kpc.shape != (LEVELS[level]**3, 3):
        raise AssertionError(f"{path}: wrong coordinate shape {gas.coords_kpc.shape}")
    for name, arr in (("coords", gas.coords_kpc), ("masses", gas.masses_code),
                      ("velocities", gas.velocities_kms), ("u", gas.internal_energy_code)):
        if np.any(~np.isfinite(arr)):
            raise AssertionError(f"{path}: non-finite {name}")
    if np.any(gas.masses_code <= 0) or np.any(gas.internal_energy_code <= 0):
        raise AssertionError(f"{path}: non-positive mass/internal energy")
    if np.min(gas.coords_kpc) <= 0 or np.max(gas.coords_kpc) >= BOX_KPC:
        raise AssertionError(f"{path}: coordinates outside open box")
    m = gas.masses_code
    com = np.sum(m[:,None] * gas.coords_kpc, axis=0) / np.sum(m)
    com_offset_pc = np.linalg.norm(com - CENTER_KPC) * 1000.0
    if com_offset_pc > 1.0e-8:
        raise AssertionError(f"{path}: gas COM offset {com_offset_pc} pc is too large")


def profile_checks(family: str, value: float, gas: GasState, level: str, turbulent_baseline: GasState | None = None) -> dict:
    rel = relative_pc(gas.coords_kpc)
    r, rhat = radial_unit_vectors(rel)
    if family == "bhl":
        cs = BHL_CS_KMS
        expected = np.zeros_like(gas.velocities_kms)
        expected[:,0] = value * cs
        err = float(np.max(np.abs(gas.velocities_kms - expected)))
        if err > 2e-12:
            raise AssertionError(f"BHL {level}: non-uniform velocity err={err}")
        mach = float(np.linalg.norm(gas.velocities_kms[0]) / cs)
        mass_msun = msun_from_code_mass(gas.masses_code)
        p = np.sum(mass_msun[:,None] * gas.velocities_kms, axis=0)
        p_expected = np.array([np.sum(mass_msun) * value * cs, 0.0, 0.0])
        if np.linalg.norm(p-p_expected) > 1e-11 * max(np.linalg.norm(p_expected), 1.0):
            raise AssertionError(f"BHL {level}: total momentum mismatch")
        return {"measured_mach": mach, "max_velocity_abs_error_kms": err}

    if family == "rotating":
        desired = rotating_velocity_field(gas.coords_kpc, value)
        mask = r < FLOW_R_FULL_PC
        dv = gas.velocities_kms - desired
        speed = np.linalg.norm(desired[mask], axis=1)
        rms_rel = math.sqrt(float(np.mean(np.sum(dv[mask]**2,axis=1)))) / max(float(np.sqrt(np.mean(speed**2))),1e-30)
        if level == "L0" and rms_rel > 1e-13:
            raise AssertionError(f"rotating {value} {level}: L0 analytic profile mismatch {rms_rel}")
        if level != "L0" and rms_rel > 0.08:
            raise AssertionError(f"rotating {value} {level}: child momentum correction too large {rms_rel}")
        mass_msun = msun_from_code_mass(gas.masses_code)
        p = np.sum(mass_msun[:,None] * gas.velocities_kms, axis=0)
        pscale = np.sum(mass_msun * np.linalg.norm(gas.velocities_kms, axis=1))
        if np.linalg.norm(p) > 2e-12 * max(pscale, 1.0):
            raise AssertionError(f"rotating {value} {level}: spurious net linear momentum")
        L = np.sum(np.cross(rel, mass_msun[:,None] * gas.velocities_kms), axis=0)
        if value == 0.0:
            Lscale = np.sum(mass_msun * r * np.linalg.norm(gas.velocities_kms,axis=1))
            if np.linalg.norm(L) > 2e-12 * max(Lscale,1.0):
                raise AssertionError(f"rotating 0 {level}: spurious angular momentum")
        else:
            if L[2] <= 0 or np.linalg.norm(L[:2]) > 2e-10 * abs(L[2]):
                raise AssertionError(f"rotating {value} {level}: angular momentum is not aligned with +z")
        return {"rms_velocity_profile_relative_correction": rms_rel}

    if family == "turbulent":
        if turbulent_baseline is None:
            raise AssertionError("missing nested spherical baseline for turbulence verification")
        dv = gas.velocities_kms - turbulent_baseline.velocities_kms
        mask = r <= TURB_NORMALIZATION_RADIUS_PC
        m = gas.masses_code[mask]
        rms = math.sqrt(float(np.sum(m*np.sum(dv[mask]**2,axis=1))/np.sum(m)))
        cs = math.sqrt(GAMMA*(GAMMA-1.0)*specific_u_from_temperature(SPHERICAL_T_K))
        mach = rms/cs
        if abs(mach-value)/value > 0.10:
            raise AssertionError(f"turbulent {value} {level}: achieved Mach {mach} differs by >10%")
        mass_msun = msun_from_code_mass(gas.masses_code)
        p = np.sum(mass_msun[:,None] * gas.velocities_kms, axis=0)
        pscale = np.sum(mass_msun * np.linalg.norm(gas.velocities_kms, axis=1))
        if np.linalg.norm(p) > 2e-11 * max(pscale, 1.0):
            raise AssertionError(f"turbulent {value} {level}: nonzero global bulk momentum")
        L = np.sum(np.cross(rel, mass_msun[:,None] * gas.velocities_kms), axis=0)
        Lscale = np.sum(mass_msun * r * np.linalg.norm(gas.velocities_kms, axis=1))
        coherence = np.linalg.norm(L) / max(Lscale, 1e-300)
        if coherence > 1.0e-3:
            raise AssertionError(f"turbulent {value} {level}: residual coherent rotation {coherence}")
        return {"measured_turbulent_mach_rle25pc": mach, "global_angular_momentum_coherence": coherence}
    raise ValueError(family)


def case_stem(family: str, tag: str, level: str) -> str:
    if family == "bhl": return f"bhl_{tag}_{level}"
    if family == "rotating": return f"rot_{tag}_{level}"
    if family == "turbulent": return f"turb_{tag}_{level}"
    raise ValueError(family)


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=HERE)
    args=p.parse_args()
    root=args.root.resolve()
    reports=[]
    total=0
    u0 = specific_u_from_temperature(SPHERICAL_T_K)
    density_fn = base_spherical_density_field
    base_velocity = base_spherical_velocity_field
    u_fn = lambda x: np.full(x.shape[0], u0, dtype=np.float64)
    turbulent_baseline = build_nested_hierarchy(density_fn, base_velocity, u_fn)
    for family, controls in EXPECTED.items():
        for tag, value in controls:
            states={}
            for level in ("L0","L1","L2"):
                stem=case_stem(family,tag,level)
                path=root/family/"ics"/level/f"{stem}.hdf5"
                if not path.exists(): raise AssertionError(f"missing {path}")
                gas,meta=load_ic(path)
                finite_positive_checks(path,gas,level)
                diag=basic_state_diagnostics(gas)
                if diag["outer_flux_shell_cells_10_12p5pc"] <= 0 or diag["shell_ffr_cells_7p5_12p5pc"] <= 0:
                    raise AssertionError(f"{path}: empty Racc diagnostic shell")
                prof=profile_checks(family,value,gas,level, turbulent_baseline[level] if family == "turbulent" else None)
                states[level]=gas
                reports.append({"family":family,"control":value,"level":level,**diag,**prof})
                total += 1
            c01=hierarchy_closure(states["L0"],states["L1"])
            c12=hierarchy_closure(states["L1"],states["L2"])
            if family == "bhl":
                density_for_split = lambda x: np.ones(x.shape[0], dtype=np.float64)
            elif family == "rotating":
                density_for_split = lambda x, f=value: rotating_density_field(x, f)
            else:
                density_for_split = base_spherical_density_field
            w01 = hierarchy_density_weight_error(states["L0"], states["L1"], density_for_split)
            w12 = hierarchy_density_weight_error(states["L1"], states["L2"], density_for_split)
            if max(w01,w12) > 5e-15:
                raise AssertionError(f"{family} {value}: child density-weight split mismatch {w01} {w12}")
            for name,c in (("L0_to_L1",c01),("L1_to_L2",c12)):
                if c["max_parent_child_mass_relative_error"] > 5e-15:
                    raise AssertionError(f"{family} {value} {name}: mass closure {c}")
                if c["max_parent_child_momentum_relative_error"] > 2e-13:
                    raise AssertionError(f"{family} {value} {name}: momentum closure {c}")
            print(f"PASS {family:9s} {value:g}: mass/momentum hierarchy {c01} {c12}")
    if total != 42:
        raise AssertionError(f"expected 42 ICs, found {total}")
    write_json(root/"verification_report.json", {"count":total,"status":"PASS","cases":reports})
    print(f"\nIC verification: PASS ({total} files)")

if __name__ == "__main__": main()
