from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import math
from typing import Callable, Dict, Iterable, Tuple

import h5py
import numpy as np

# -----------------------------
# Physical/code-unit constants
# -----------------------------
G_CGS = 6.67430e-8
K_B_CGS = 1.380649e-16
M_P_CGS = 1.67262192369e-24
M_SUN_CGS = 1.98847e33
PC_CGS = 3.0856775814913673e18
KPC_CGS = 1.0e3 * PC_CGS
YR_CGS = 365.25 * 86400.0
KM_CGS = 1.0e5

UNIT_LENGTH_CGS = KPC_CGS
UNIT_MASS_CGS = 1.0e10 * M_SUN_CGS
UNIT_VELOCITY_CGS = KM_CGS
UNIT_ENERGY_PER_MASS_CGS = UNIT_VELOCITY_CGS**2

BOX_PC = 100.0
BOX_KPC = BOX_PC / 1000.0
CENTER_KPC = np.array([0.5 * BOX_KPC] * 3, dtype=np.float64)
BH_MASS_MSUN = 1.0e5
BH_MASS_CODE = BH_MASS_MSUN / 1.0e10
BH_ID = np.uint64(1000000001)
R_ACC_PC = 12.5
R_FEEDBACK_PC = 25.0
R_CORE_PC = 0.5
FLOW_R_FULL_PC = 35.0
FLOW_R_TAPER_PC = 40.0
MDOT_TRUE_MSUN_PER_YR = 1.0e-4
MU = 1.22
GAMMA = 5.0 / 3.0
X_H = 0.76
SPHERICAL_T_K = 100.0
BHL_CS_KMS = 5.0
BHL_FIDUCIAL_MACH = 1.0
TURB_SEED = 20261006
TURB_KMAX = 4
TURB_NORMALIZATION_RADIUS_PC = R_FEEDBACK_PC

LEVELS = {"L0": 16, "L1": 32, "L2": 64}


@dataclass
class GasState:
    coords_kpc: np.ndarray
    masses_code: np.ndarray
    velocities_kms: np.ndarray
    internal_energy_code: np.ndarray

    def copy(self) -> "GasState":
        return GasState(
            self.coords_kpc.copy(),
            self.masses_code.copy(),
            self.velocities_kms.copy(),
            self.internal_energy_code.copy(),
        )


def code_mass_from_msun(mass_msun: float | np.ndarray) -> float | np.ndarray:
    return np.asarray(mass_msun) / 1.0e10


def msun_from_code_mass(mass_code: float | np.ndarray) -> float | np.ndarray:
    return np.asarray(mass_code) * 1.0e10


def specific_u_from_temperature(T_K: float, mu: float = MU, gamma: float = GAMMA) -> float:
    u_cgs = K_B_CGS * T_K / ((gamma - 1.0) * mu * M_P_CGS)
    return u_cgs / UNIT_ENERGY_PER_MASS_CGS


def temperature_from_sound_speed(cs_kms: float, mu: float = MU, gamma: float = GAMMA) -> float:
    cs_cgs = cs_kms * KM_CGS
    return mu * M_P_CGS * cs_cgs**2 / (gamma * K_B_CGS)


def specific_u_from_sound_speed(cs_kms: float, gamma: float = GAMMA) -> float:
    return cs_kms**2 / (gamma * (gamma - 1.0))


def grav_g_kpc_kms2_per_msun() -> float:
    # G in kpc (km/s)^2 / Msun
    return G_CGS * M_SUN_CGS / KPC_CGS / (KM_CGS**2)


G_KPC_KMS2_PER_MSUN = grav_g_kpc_kms2_per_msun()
G_PC_KMS2_PER_MSUN = 1000.0 * G_KPC_KMS2_PER_MSUN


def regular_l0_positions(nside: int = 16) -> np.ndarray:
    dx = BOX_KPC / nside
    q = (np.arange(nside, dtype=np.float64) + 0.5) * dx
    x, y, z = np.meshgrid(q, q, q, indexing="ij")
    return np.column_stack([x.ravel(), y.ravel(), z.ravel()])


def relative_pc(coords_kpc: np.ndarray) -> np.ndarray:
    return (coords_kpc - CENTER_KPC[None, :]) * 1000.0


def smoothstep01(x: np.ndarray) -> np.ndarray:
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


def flow_taper(r_pc: np.ndarray) -> np.ndarray:
    out = np.ones_like(r_pc, dtype=np.float64)
    mid = (r_pc > FLOW_R_FULL_PC) & (r_pc < FLOW_R_TAPER_PC)
    out[r_pc >= FLOW_R_TAPER_PC] = 0.0
    if np.any(mid):
        x = (r_pc[mid] - FLOW_R_FULL_PC) / (FLOW_R_TAPER_PC - FLOW_R_FULL_PC)
        out[mid] = 1.0 - smoothstep01(x)
    return out


def spherical_profile_density_cgs(r_pc: np.ndarray, radial_speed_kms: np.ndarray) -> np.ndarray:
    """Constant-Mdot density inside 35 pc, tapered to a 1%-density ambient by 40 pc."""
    r_eff = np.maximum(r_pc, R_CORE_PC)
    v_eff_cgs = np.maximum(np.abs(radial_speed_kms), 1.0e-30) * KM_CGS
    mdot_cgs = MDOT_TRUE_MSUN_PER_YR * M_SUN_CGS / YR_CGS
    rho_flow = mdot_cgs / (4.0 * np.pi * (r_eff * PC_CGS) ** 2 * v_eff_cgs)

    # Define the ambient level from the un-tapered solution at the outer-flow radius.
    # radial_speed_kms is model-dependent, so infer the local r^{-1/2} normalization
    # from each point and evaluate it at FLOW_R_FULL_PC.
    v_outer = np.abs(radial_speed_kms) * np.sqrt(np.maximum(r_eff, 1e-30) / FLOW_R_FULL_PC)
    rho_outer = mdot_cgs / (
        4.0 * np.pi * (FLOW_R_FULL_PC * PC_CGS) ** 2 * np.maximum(v_outer, 1e-30) * KM_CGS
    )
    rho_ambient = 0.01 * rho_outer

    taper = flow_taper(r_pc)
    return taper * rho_flow + (1.0 - taper) * rho_ambient


def radial_unit_vectors(rel_pc: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    r = np.linalg.norm(rel_pc, axis=1)
    rhat = np.zeros_like(rel_pc)
    good = r > 0
    rhat[good] = rel_pc[good] / r[good, None]
    return r, rhat


def azimuthal_unit_vectors(rel_pc: np.ndarray) -> np.ndarray:
    x = rel_pc[:, 0]
    y = rel_pc[:, 1]
    R = np.hypot(x, y)
    phihat = np.zeros_like(rel_pc)
    good = R > 0
    phihat[good, 0] = -y[good] / R[good]
    phihat[good, 1] = x[good] / R[good]
    return phihat


def rotating_radial_speed_kms(r_pc: np.ndarray, frot: float) -> np.ndarray:
    r_eff = np.maximum(r_pc, R_CORE_PC)
    return -np.sqrt((2.0 - frot * frot) * G_PC_KMS2_PER_MSUN * BH_MASS_MSUN / r_eff)


def rotating_velocity_field(coords_kpc: np.ndarray, frot: float) -> np.ndarray:
    rel = relative_pc(coords_kpc)
    r, rhat = radial_unit_vectors(rel)
    r_eff = np.maximum(r, R_CORE_PC)
    phihat = azimuthal_unit_vectors(rel)
    vr = rotating_radial_speed_kms(r, frot)
    vphi = frot * np.sqrt(G_PC_KMS2_PER_MSUN * BH_MASS_MSUN / r_eff)
    taper = flow_taper(r)
    return taper[:, None] * (vr[:, None] * rhat + vphi[:, None] * phihat)


def rotating_density_field(coords_kpc: np.ndarray, frot: float) -> np.ndarray:
    rel = relative_pc(coords_kpc)
    r = np.linalg.norm(rel, axis=1)
    vr = rotating_radial_speed_kms(r, frot)
    return spherical_profile_density_cgs(r, vr)


def base_spherical_density_field(coords_kpc: np.ndarray) -> np.ndarray:
    return rotating_density_field(coords_kpc, 0.0)


def base_spherical_velocity_field(coords_kpc: np.ndarray) -> np.ndarray:
    return rotating_velocity_field(coords_kpc, 0.0)


def masses_from_density_on_l0(coords_kpc: np.ndarray, density_cgs: np.ndarray, nside: int = 16) -> np.ndarray:
    dx_cgs = (BOX_KPC / nside) * KPC_CGS
    mass_g = density_cgs * dx_cgs**3
    return mass_g / UNIT_MASS_CGS


def split_state(
    parent: GasState,
    parent_dx_kpc: float,
    density_fn: Callable[[np.ndarray], np.ndarray],
    velocity_fn: Callable[[np.ndarray], np.ndarray],
    internal_energy_fn: Callable[[np.ndarray], np.ndarray],
) -> GasState:
    """Split every parent into eight children while conserving parent mass and momentum exactly."""
    signs = np.array(
        [[sx, sy, sz] for sx in (-1.0, 1.0) for sy in (-1.0, 1.0) for sz in (-1.0, 1.0)],
        dtype=np.float64,
    )
    offsets = signs * (parent_dx_kpc / 4.0)
    npar = parent.coords_kpc.shape[0]
    child_coords = (parent.coords_kpc[:, None, :] + offsets[None, :, :]).reshape(-1, 3)

    rho = np.asarray(density_fn(child_coords), dtype=np.float64).reshape(npar, 8)
    if np.any(~np.isfinite(rho)) or np.any(rho <= 0):
        raise ValueError("non-positive/non-finite child target density")
    weights = rho / np.sum(rho, axis=1)[:, None]
    child_mass = (parent.masses_code[:, None] * weights).reshape(-1)

    desired_vel = np.asarray(velocity_fn(child_coords), dtype=np.float64).reshape(npar, 8, 3)
    mass8 = child_mass.reshape(npar, 8)
    vmean = np.sum(mass8[:, :, None] * desired_vel, axis=1) / parent.masses_code[:, None]
    correction = parent.velocities_kms - vmean
    child_vel = (desired_vel + correction[:, None, :]).reshape(-1, 3)

    child_u = np.asarray(internal_energy_fn(child_coords), dtype=np.float64).reshape(-1)
    if np.any(~np.isfinite(child_u)) or np.any(child_u <= 0):
        raise ValueError("non-positive/non-finite child internal energy")

    return GasState(child_coords, child_mass, child_vel, child_u)


def build_nested_hierarchy(
    density_fn: Callable[[np.ndarray], np.ndarray],
    velocity_fn: Callable[[np.ndarray], np.ndarray],
    internal_energy_fn: Callable[[np.ndarray], np.ndarray],
) -> Dict[str, GasState]:
    l0_coords = regular_l0_positions(16)
    l0_rho = density_fn(l0_coords)
    l0_mass = masses_from_density_on_l0(l0_coords, l0_rho, 16)
    l0_vel = velocity_fn(l0_coords)
    l0_u = internal_energy_fn(l0_coords)
    l0 = GasState(l0_coords, l0_mass, l0_vel, l0_u)
    l1 = split_state(l0, BOX_KPC / 16.0, density_fn, velocity_fn, internal_energy_fn)
    l2 = split_state(l1, BOX_KPC / 32.0, density_fn, velocity_fn, internal_energy_fn)
    return {"L0": l0, "L1": l1, "L2": l2}


def write_arepo_ic(path: Path, gas: GasState, metadata: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ngas = gas.coords_kpc.shape[0]
    ids = np.arange(1, ngas + 1, dtype=np.uint64)
    n_this = np.array([ngas, 0, 0, 0, 0, 1], dtype=np.uint32)

    with h5py.File(path, "w") as f:
        h = f.create_group("Header")
        h.attrs["NumPart_ThisFile"] = n_this
        h.attrs["NumPart_Total"] = n_this
        h.attrs["NumPart_Total_HighWord"] = np.zeros(6, dtype=np.uint32)
        h.attrs["MassTable"] = np.zeros(6, dtype=np.float64)
        h.attrs["Time"] = 0.0
        h.attrs["Redshift"] = 0.0
        h.attrs["BoxSize"] = BOX_KPC
        h.attrs["NumFilesPerSnapshot"] = 1
        h.attrs["Omega0"] = 0.0
        h.attrs["OmegaLambda"] = 0.0
        h.attrs["HubbleParam"] = 1.0
        h.attrs["Flag_Sfr"] = 0
        h.attrs["Flag_Cooling"] = 0
        h.attrs["Flag_StellarAge"] = 0
        h.attrs["Flag_Metals"] = 0
        h.attrs["Flag_Feedback"] = 0
        h.attrs["Flag_DoublePrecision"] = 1

        g = f.create_group("PartType0")
        g.create_dataset("Coordinates", data=np.asarray(gas.coords_kpc, dtype=np.float64))
        g.create_dataset("Velocities", data=np.asarray(gas.velocities_kms, dtype=np.float64))
        g.create_dataset("ParticleIDs", data=ids)
        g.create_dataset("Masses", data=np.asarray(gas.masses_code, dtype=np.float64))
        g.create_dataset("InternalEnergy", data=np.asarray(gas.internal_energy_code, dtype=np.float64))

        bh = f.create_group("PartType5")
        bh.create_dataset("Coordinates", data=CENTER_KPC.reshape(1, 3))
        bh.create_dataset("Velocities", data=np.zeros((1, 3), dtype=np.float64))
        bh.create_dataset("ParticleIDs", data=np.array([BH_ID], dtype=np.uint64))
        bh.create_dataset("Masses", data=np.array([BH_MASS_CODE], dtype=np.float64))

        m = f.create_group("ICMetadata")
        m.attrs["Generator"] = "bh_accretion_flow_tests/generate_all_ics.py"
        m.attrs["SchemaVersion"] = 1
        m.attrs["Family"] = metadata["family"]
        m.attrs["Level"] = metadata["level"]
        m.attrs["ControlName"] = metadata["control_name"]
        m.attrs["ControlValue"] = metadata["control_value"]


def write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True)
        f.write("\n")


def cell_dx_pc(level: str) -> float:
    return BOX_PC / LEVELS[level]


def basic_state_diagnostics(gas: GasState) -> dict:
    rel_pc = relative_pc(gas.coords_kpc)
    r_pc = np.linalg.norm(rel_pc, axis=1)
    mass_msun = msun_from_code_mass(gas.masses_code)
    total = float(np.sum(mass_msun))
    com = np.sum(mass_msun[:, None] * gas.coords_kpc, axis=0) / total
    p = np.sum(mass_msun[:, None] * gas.velocities_kms, axis=0)
    L = np.sum(np.cross(rel_pc, mass_msun[:, None] * gas.velocities_kms), axis=0)
    shell_outer = (r_pc >= 10.0) & (r_pc < R_ACC_PC)
    shell_ffr = (r_pc >= 0.6 * R_ACC_PC) & (r_pc < R_ACC_PC)
    return {
        "ngas": int(gas.coords_kpc.shape[0]),
        "total_gas_mass_msun": total,
        "center_of_mass_kpc": com.tolist(),
        "net_momentum_msun_kms": p.tolist(),
        "net_angular_momentum_msun_pc_kms": L.tolist(),
        "outer_flux_shell_cells_10_12p5pc": int(np.count_nonzero(shell_outer)),
        "shell_ffr_cells_7p5_12p5pc": int(np.count_nonzero(shell_ffr)),
    }


def shell_flux_msun_per_yr(gas: GasState, rlo_pc: float, rhi_pc: float) -> Tuple[float, float, int]:
    rel = relative_pc(gas.coords_kpc)
    r, rhat = radial_unit_vectors(rel)
    mask = (r >= rlo_pc) & (r < rhi_pc)
    if not np.any(mask):
        return float("nan"), float("nan"), 0
    vr_kms = np.sum(gas.velocities_kms * rhat, axis=1)
    mass_msun = msun_from_code_mass(gas.masses_code)
    # -sum(m v_r)/dr, converting km/s/pc -> 1/yr.
    conv = KM_CGS / PC_CGS * YR_CGS
    net = -np.sum(mass_msun[mask] * vr_kms[mask]) / (rhi_pc - rlo_pc) * conv
    inward = np.sum(mass_msun[mask] * np.maximum(-vr_kms[mask], 0.0)) / (rhi_pc - rlo_pc) * conv
    return float(net), float(inward), int(np.count_nonzero(mask))


def vphi_proxy_kms(gas: GasState, radius_pc: float = R_ACC_PC) -> float:
    rel = relative_pc(gas.coords_kpc)
    r = np.linalg.norm(rel, axis=1)
    mask = r < radius_pc
    if not np.any(mask):
        return float("nan")
    m = gas.masses_code[mask]
    Lvec = np.sum(np.cross(rel[mask], m[:, None] * gas.velocities_kms[mask]), axis=0)
    return float(np.linalg.norm(Lvec) / (np.sum(m) * radius_pc))


class SolenoidalTurbulence:
    def __init__(self, seed: int = TURB_SEED, kmax: int = TURB_KMAX):
        self.seed = int(seed)
        self.kmax = int(kmax)
        rng = np.random.default_rng(self.seed)
        modes = []
        for nx in range(-kmax, kmax + 1):
            for ny in range(-kmax, kmax + 1):
                for nz in range(-kmax, kmax + 1):
                    if nx == ny == nz == 0:
                        continue
                    n = np.array([nx, ny, nz], dtype=np.float64)
                    kmag = np.linalg.norm(n)
                    if kmag < 1.0 or kmag > kmax + 1e-12:
                        continue
                    # Keep one representative of each +/- pair.
                    first = next(int(v) for v in (nx, ny, nz) if v != 0)
                    if first < 0:
                        continue
                    a = rng.normal(size=3)
                    a -= np.dot(a, n) / np.dot(n, n) * n
                    an = np.linalg.norm(a)
                    if an < 1e-14:
                        continue
                    a /= an
                    amp = kmag ** -2.0  # sqrt(P_v) for P_v ~ k^-4
                    phase = rng.uniform(0.0, 2.0 * np.pi)
                    modes.append((np.array([nx, ny, nz], dtype=np.int32), amp * a, phase))
        self.modes = modes
        self.bulk = np.zeros(3)
        self.omega = np.zeros(3)
        self.scale = 1.0

    def raw(self, coords_kpc: np.ndarray) -> np.ndarray:
        x = coords_kpc / BOX_KPC
        v = np.zeros((coords_kpc.shape[0], 3), dtype=np.float64)
        for n, avec, phase in self.modes:
            arg = 2.0 * np.pi * (x @ n.astype(np.float64)) + phase
            v += np.cos(arg)[:, None] * avec[None, :]
        return v

    def corrected_unscaled(self, coords_kpc: np.ndarray) -> np.ndarray:
        rel = relative_pc(coords_kpc)
        return self.raw(coords_kpc) - self.bulk[None, :] - np.cross(self.omega[None, :], rel)

    def __call__(self, coords_kpc: np.ndarray) -> np.ndarray:
        return self.scale * self.corrected_unscaled(coords_kpc)

    def calibrate(self, coords_kpc: np.ndarray, masses_code: np.ndarray, target_rms_kms: float,
                  radius_pc: float = TURB_NORMALIZATION_RADIUS_PC) -> None:
        rel = relative_pc(coords_kpc)
        r = np.linalg.norm(rel, axis=1)
        mask_rms = r <= radius_pc
        if np.count_nonzero(mask_rms) < 8:
            raise ValueError("too few cells in turbulence normalization region")

        # Remove the *global* mass-weighted bulk velocity and rigid-body rotation,
        # so the turbulent perturbation adds neither net linear momentum nor
        # coherent global angular momentum to the spherical inflow.  The target
        # RMS amplitude is then imposed only in the documented normalization region.
        m_all = masses_code
        vv_all = self.raw(coords_kpc)
        mtot_all = np.sum(m_all)
        bulk = np.sum(m_all[:, None] * vv_all, axis=0) / mtot_all
        vv0 = vv_all - bulk[None, :]
        inertia = np.zeros((3, 3), dtype=np.float64)
        for mi, ri in zip(m_all, rel):
            r2 = float(np.dot(ri, ri))
            inertia += mi * (r2 * np.eye(3) - np.outer(ri, ri))
        L = np.sum(np.cross(rel, m_all[:, None] * vv0), axis=0)
        omega = np.linalg.solve(inertia, L)
        vv1 = vv0 - np.cross(omega[None, :], rel)
        bulk2 = np.sum(m_all[:, None] * vv1, axis=0) / mtot_all

        self.bulk = bulk + bulk2
        self.omega = omega
        vv_corr = self.corrected_unscaled(coords_kpc[mask_rms])
        m_rms = masses_code[mask_rms]
        rms = math.sqrt(float(np.sum(m_rms * np.sum(vv_corr**2, axis=1)) / np.sum(m_rms)))
        if not (rms > 0):
            raise ValueError("zero turbulence RMS")
        self.scale = target_rms_kms / rms

    def mode_table(self) -> list:
        out = []
        for n, a, ph in self.modes:
            out.append({
                "n": [int(x) for x in n],
                "amplitude_vector_unscaled": [float(x) for x in a],
                "phase_rad": float(ph),
            })
        return out


def turbulence_rms(gas: GasState, turb_fn: SolenoidalTurbulence, radius_pc: float) -> float:
    rel = relative_pc(gas.coords_kpc)
    r = np.linalg.norm(rel, axis=1)
    mask = r <= radius_pc
    m = gas.masses_code[mask]
    vt = turb_fn(gas.coords_kpc[mask])
    # The actual hierarchy velocities contain parent-momentum corrections, so this
    # helper is only for the underlying continuous realization.
    return math.sqrt(float(np.sum(m * np.sum(vt**2, axis=1)) / np.sum(m)))


def actual_turbulent_component_rms(gas: GasState, radius_pc: float) -> float:
    rel = relative_pc(gas.coords_kpc)
    r = np.linalg.norm(rel, axis=1)
    mask = r <= radius_pc
    base = base_spherical_velocity_field(gas.coords_kpc)
    dv = gas.velocities_kms - base
    m = gas.masses_code
    return math.sqrt(float(np.sum(m[mask] * np.sum(dv[mask]**2, axis=1)) / np.sum(m[mask])))


def hierarchy_closure(parent: GasState, child: GasState) -> dict:
    npar = parent.coords_kpc.shape[0]
    if child.coords_kpc.shape[0] != 8 * npar:
        raise ValueError("child count is not 8x parent count")
    cm = child.masses_code.reshape(npar, 8)
    cv = child.velocities_kms.reshape(npar, 8, 3)
    msum = np.sum(cm, axis=1)
    psum = np.sum(cm[:, :, None] * cv, axis=1)
    ppar = parent.masses_code[:, None] * parent.velocities_kms
    mass_rel = np.max(np.abs(msum - parent.masses_code) / np.maximum(np.abs(parent.masses_code), 1e-300))
    pscale = np.maximum(np.linalg.norm(ppar, axis=1), np.sum(cm * np.linalg.norm(cv, axis=2), axis=1))
    pdiff = np.linalg.norm(psum - ppar, axis=1)
    mom_rel = np.max(pdiff / np.maximum(pscale, 1e-300))
    return {"max_parent_child_mass_relative_error": float(mass_rel),
            "max_parent_child_momentum_relative_error": float(mom_rel)}


def hierarchy_density_weight_error(parent: GasState, child: GasState, density_fn: Callable[[np.ndarray], np.ndarray]) -> float:
    npar = parent.coords_kpc.shape[0]
    if child.coords_kpc.shape[0] != 8 * npar:
        raise ValueError("child count is not 8x parent count")
    cm = child.masses_code.reshape(npar, 8)
    rho = np.asarray(density_fn(child.coords_kpc), dtype=np.float64).reshape(npar, 8)
    expected = parent.masses_code[:, None] * rho / np.sum(rho, axis=1)[:, None]
    return float(np.max(np.abs(cm - expected) / np.maximum(np.abs(expected), 1e-300)))
