#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
G = 4.30091e-6  # kpc (km/s)^2 Msun^-1
GAMMA = 5.0 / 3.0
BOX_KPC = 200.0
CENTER = np.array([BOX_KPC / 2.0] * 3)
SEED = 20261006

# Compact isolated-galaxy benchmark.  This is intentionally not a precision
# equilibrium galaxy model; the BH remains off for 20 Myr so transient IC
# relaxation occurs before any gas is captured.
M_HALO = 5.0e10
A_HALO = 15.0
N_HALO = 80_000

M_STELLAR_DISK = 3.0e9
R_STELLAR_DISK = 2.2
Z_STELLAR_DISK = 0.22
N_STELLAR_DISK = 30_000

M_BULGE = 5.0e8
A_BULGE = 0.55
N_BULGE = 10_000

M_GAS_EXT = 6.5e8
R_GAS_EXT = 2.8
Z_GAS_EXT = 0.16
N_GAS_EXT = 40_000

M_GAS_NUC = 1.5e8
R_GAS_NUC = 0.35
Z_GAS_NUC = 0.05
N_GAS_NUC = 20_000

M_BH = 1.0e6

SETTLE_MYR = 20.0
END_MYR = 80.0
RATE_CODE_TO_MSUN_YR = 10.2271202634
MYR_PER_CODE = 977.792354298


def sample_exponential_disk(rng, n, rd, zd, rmax):
    out_r = []
    while sum(len(x) for x in out_r) < n:
        r = rng.gamma(shape=2.0, scale=rd, size=max(n, 4096))
        r = r[r < rmax]
        out_r.append(r)
    R = np.concatenate(out_r)[:n]
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    z = rng.exponential(zd, n) * rng.choice([-1.0, 1.0], n)
    xyz = np.column_stack((R * np.cos(phi), R * np.sin(phi), z))
    return xyz, R, phi


def sample_hernquist(rng, n, a, rmax):
    chunks = []
    while sum(len(x) for x in chunks) < n:
        u = rng.uniform(1.0e-9, 1.0 - 1.0e-9, max(n, 4096))
        y = np.sqrt(u)
        r = a * y / (1.0 - y)
        chunks.append(r[r < rmax])
    r = np.concatenate(chunks)[:n]
    cost = rng.uniform(-1.0, 1.0, n)
    sint = np.sqrt(1.0 - cost * cost)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    xyz = np.column_stack((r * sint * np.cos(phi),
                           r * sint * np.sin(phi),
                           r * cost))
    return xyz, r


def disk_enclosed(m, R, rd):
    x = np.maximum(R / rd, 0.0)
    return m * (1.0 - (1.0 + x) * np.exp(-x))


def hernquist_enclosed(m, r, a):
    return m * r * r / (r + a) ** 2


def vcirc(R):
    R = np.maximum(np.asarray(R, dtype=np.float64), 0.03)
    menc = (
        hernquist_enclosed(M_HALO, R, A_HALO)
        + hernquist_enclosed(M_BULGE, R, A_BULGE)
        + disk_enclosed(M_STELLAR_DISK, R, R_STELLAR_DISK)
        + disk_enclosed(M_GAS_EXT, R, R_GAS_EXT)
        + disk_enclosed(M_GAS_NUC, R, R_GAS_NUC)
        + M_BH
    )
    return np.sqrt(G * menc / R)


def cylindrical_velocities(rng, R, phi, rotation_fraction,
                           sigma_r, sigma_phi, sigma_z, radial_mean=0.0):
    v0 = vcirc(R)
    vr = radial_mean + rng.normal(0.0, sigma_r, len(R))
    vp = rotation_fraction * v0 + rng.normal(0.0, sigma_phi, len(R))
    vz = rng.normal(0.0, sigma_z, len(R))
    vx = vr * np.cos(phi) - vp * np.sin(phi)
    vy = vr * np.sin(phi) + vp * np.cos(phi)
    return np.column_stack((vx, vy, vz))


def isotropic_velocities(rng, xyz, scale=0.55):
    r = np.linalg.norm(xyz, axis=1)
    sigma = scale * vcirc(np.maximum(r, 0.05))
    vel = rng.normal(size=(len(r), 3)) * sigma[:, None]
    return vel


def write_group(f, ptype, pos, vel, mass_msun, ids, u=None):
    g = f.create_group(f"PartType{ptype}")
    g.create_dataset("Coordinates", data=np.asarray(pos, dtype=np.float64))
    g.create_dataset("Velocities", data=np.asarray(vel, dtype=np.float64))
    g.create_dataset("ParticleIDs", data=np.asarray(ids, dtype=np.uint64))
    mass = np.asarray(mass_msun, dtype=np.float64)
    if mass.ndim == 0:
        mass = np.full(len(pos), float(mass), dtype=np.float64)
    g.create_dataset("Masses", data=mass / 1.0e10)
    if u is not None:
        g.create_dataset("InternalEnergy", data=np.asarray(u, dtype=np.float64))


def main():
    rng = np.random.default_rng(SEED)
    icdir = HERE / "ics"
    outdir = HERE / "outputs" / "galaxy_env_shell"
    icdir.mkdir(parents=True, exist_ok=True)
    (HERE / "outputlists").mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    # Gas: extended + deliberately well-resolved nuclear component.
    xge, Rge, phige = sample_exponential_disk(rng, N_GAS_EXT, R_GAS_EXT, Z_GAS_EXT, 15.0)
    xgn, Rgn, phign = sample_exponential_disk(rng, N_GAS_NUC, R_GAS_NUC, Z_GAS_NUC, 2.0)
    vge = cylindrical_velocities(rng, Rge, phige, 0.93, 6.0, 6.0, 4.0, radial_mean=-2.0)
    vgn = cylindrical_velocities(rng, Rgn, phign, 0.85, 12.0, 10.0, 8.0, radial_mean=-6.0)
    gas_pos = np.vstack((xge, xgn)) + CENTER
    gas_vel = np.vstack((vge, vgn))
    gas_mass = np.concatenate((
        np.full(N_GAS_EXT, M_GAS_EXT / N_GAS_EXT),
        np.full(N_GAS_NUC, M_GAS_NUC / N_GAS_NUC),
    ))
    cs = np.concatenate((np.full(N_GAS_EXT, 12.0), np.full(N_GAS_NUC, 20.0)))
    gas_u = cs * cs / (GAMMA * (GAMMA - 1.0))

    # Stellar disk and bulge.
    xsd, Rsd, phisd = sample_exponential_disk(
        rng, N_STELLAR_DISK, R_STELLAR_DISK, Z_STELLAR_DISK, 12.0
    )
    vsd = cylindrical_velocities(rng, Rsd, phisd, 0.97, 22.0, 16.0, 12.0)
    xb, _ = sample_hernquist(rng, N_BULGE, A_BULGE, 8.0)
    vb = isotropic_velocities(rng, xb, scale=0.55)
    star_pos = np.vstack((xsd, xb)) + CENTER
    star_vel = np.vstack((vsd, vb))
    star_mass = np.concatenate((
        np.full(N_STELLAR_DISK, M_STELLAR_DISK / N_STELLAR_DISK),
        np.full(N_BULGE, M_BULGE / N_BULGE),
    ))

    # Live Hernquist-like dark halo.
    xh, _ = sample_hernquist(rng, N_HALO, A_HALO, 80.0)
    vh = isotropic_velocities(rng, xh, scale=0.60)
    halo_pos = xh + CENTER
    halo_vel = vh
    halo_mass = np.full(N_HALO, M_HALO / N_HALO)

    # Remove small net COM drifts from each collisionless component.
    halo_vel -= np.average(halo_vel, axis=0, weights=halo_mass)
    star_vel -= np.average(star_vel, axis=0, weights=star_mass)

    n0 = len(gas_pos)
    n1 = len(halo_pos)
    n4 = len(star_pos)
    n5 = 1
    counts = np.array([n0, n1, 0, 0, n4, n5], dtype=np.uint32)

    ids0 = np.arange(1, n0 + 1, dtype=np.uint64)
    ids1 = np.arange(ids0[-1] + 1, ids0[-1] + 1 + n1, dtype=np.uint64)
    ids4 = np.arange(ids1[-1] + 1, ids1[-1] + 1 + n4, dtype=np.uint64)
    ids5 = np.array([1_000_000_001], dtype=np.uint64)

    path = icdir / "isolated_disk_galaxy.hdf5"
    with h5py.File(path, "w") as f:
        h = f.create_group("Header")
        h.attrs["NumPart_ThisFile"] = counts
        h.attrs["NumPart_Total"] = counts
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
        h.attrs["Flag_Entropy_ICs"] = 0

        write_group(f, 0, gas_pos, gas_vel, gas_mass, ids0, u=gas_u)
        write_group(f, 1, halo_pos, halo_vel, halo_mass, ids1)
        write_group(f, 4, star_pos, star_vel, star_mass, ids4)
        write_group(f, 5, CENTER[None, :], np.zeros((1, 3)), M_BH, ids5)

    snap_myr = [5, 10, 15, 19, 20, 21, 25, 30, 40, 50, 60, 70, 80]
    with (HERE / "outputlists" / "galaxy_times.txt").open("w") as f:
        for t in snap_myr:
            f.write(f"{t / MYR_PER_CODE:.17g} 2\n")

    tmax = END_MYR / MYR_PER_CODE
    maxdt = 0.25 / MYR_PER_CODE
    statdt = 1.0 / MYR_PER_CODE

    param = f"""% Isolated live-galaxy accretion comparison.
% BH gravity is present from t=0. Gas capture is disabled until {SETTLE_MYR:g} Myr.
InitCondFile                            ics/isolated_disk_galaxy
OutputDir                               outputs/galaxy_env_shell
SnapshotFileBase                        snap
OutputListFilename                      outputlists/galaxy_times.txt
ICFormat                                3
SnapFormat                              3

TimeLimitCPU                            172800
CpuTimeBetRestartFile                   1.0e30
FlushCpuTimeDiff                        120
ResubmitOn                              0
ResubmitCommand                         none
MaxMemSize                              3000

TimeBegin                               0.0
TimeMax                                 {tmax:.17g}
ComovingIntegrationOn                   0
PeriodicBoundariesOn                    1
CoolingOn                               0
StarformationOn                         0
Omega0                                  0.0
OmegaLambda                             0.0
OmegaBaryon                             0.0
HubbleParam                             1.0
BoxSize                                 {BOX_KPC:.17g}

OutputListOn                            1
TimeBetSnapshot                         1.0
TimeOfFirstSnapshot                     1.0
TimeBetStatistics                       {statdt:.17g}
NumFilesPerSnapshot                     1
NumFilesWrittenInParallel               1

TypeOfTimestepCriterion                 0
ErrTolIntAccuracy                       0.012
CourantFac                              0.3
MaxSizeTimestep                         {maxdt:.17g}
MinSizeTimestep                         1.0e-12

InitGasTemp                             10000.0
MinGasTemp                              100.0
MinimumDensityOnStartUp                 1.0e-20
LimitUBelowThisDensity                  0.0
LimitUBelowCertainDensityToThisValue    0.0
MinEgySpec                              0.0

TypeOfOpeningCriterion                  1
ErrTolTheta                             0.7
ErrTolForceAcc                          0.0025
MultipleDomains                         4
TopNodeFactor                           2.5
ActivePartFracForNewDomainDecomp        0.01
DesNumNgb                               64
MaxNumNgbDeviation                      4

UnitLength_in_cm                        3.085678e21
UnitMass_in_g                           1.989e43
UnitVelocity_in_cm_per_s                1.0e5
GravityConstantInternal                 0

SofteningComovingType0                  0.01
SofteningComovingType1                  0.05
SofteningMaxPhysType0                   0.01
SofteningMaxPhysType1                   0.05
GasSoftFactor                           2.5
SofteningTypeOfPartType0                0
SofteningTypeOfPartType1                1
SofteningTypeOfPartType2                1
SofteningTypeOfPartType3                1
SofteningTypeOfPartType4                1
SofteningTypeOfPartType5                1
MinimumComovingHydroSoftening           0.0025
AdaptiveHydroSofteningSpacing           1.2

CellShapingSpeed                        0.5
CellShapingFactor                       1.0

BHAccretionRadius                       0.30
BHFeedbackRadius                        0.60
BHSeedHaloMassMsun                      2.0e12
BHSeedMassMsun                          1.0e6
BHSeedMinRedshift                       100.0
BHSeedMaxDonorFraction                  0.5

BHFreeFallA                             1.0e-3
BHFreeFallAlpha                         0.0
BHUseFullCentralMass                    1
BHMaxSinkFraction                       0.10

BHDiskTimeMyr                           42.0
BHDiskTimeExponent                      0.4
BHRHotMaxInRs                           1200.0
BHRadiativeEfficiency                   0.1
BHWindBurstFactor                       0.01
BHJetBurstFactor                        0.10
BHUseCentralBindingTerm                 1
BHWindConeAngleDeg                      45.0
BHJetConeAngleDeg                       15.0
BHMinCoherence                          0.30
BHJetDirectionTimeFactor                1.0
BHMinTargetsPerLobe                     8
BHMinActiveTargetMassFrac               0.5
BHMaxPacketsPerStep                     4
BHDMNeighbours                          32
BHInternalTimestepFactor                1.0e30

% Selected physical path: corrected shell supply -> unresolved reservoir.
BHBenchmarkAccretionModel               5
BHBenchmarkAccretionTarget              0
BHBenchmarkBoostMode                    1
BHBenchmarkBoostAlpha                   100.0
BHBenchmarkBoostDensityThreshold        0.1
BHBenchmarkBoostBeta                    2.0
BHBenchmarkAMViscosity                  6.283185307179586
BHBenchmarkRadiativeEfficiency          0.2
BHBenchmarkEddingtonFactor              10.0
BHBenchmarkEnvThermalFloor              0.4
BHBenchmarkEnvRotationBeta              1.0
BHBenchmarkStartTimeMyr                 {SETTLE_MYR:.17g}

% Feedback intentionally disabled: compare accretion prescriptions first.
BHBenchmarkFeedbackModel                0
BHBenchmarkTNGChi0                      0.002
BHBenchmarkTNGChiBeta                   2.0
BHBenchmarkTNGChiMax                    0.1
BHBenchmarkTNGThermalCoupling           0.1
BHBenchmarkTNGKineticMaxEfficiency      0.2
BHBenchmarkTNGKineticDensityFactor      0.05
BHBenchmarkTNGSFThresholdNH             0.1065
BHBenchmarkTNGKineticBurstFactor        20.0
"""
    (HERE / "param_galaxy.txt").write_text(param)

    def set_param(text: str, key: str, value: str) -> str:
        lines = text.splitlines()
        found = False
        for i, line in enumerate(lines):
            if line.startswith(key):
                lines[i] = f"{key:<40s}{value}"
                found = True
                break
        if not found:
            raise KeyError(key)
        return "\n".join(lines) + "\n"

    # Stage 1: one common BH-off settling run.
    settle_times = [5, 10, 15, 20]
    with (HERE / "outputlists" / "settle_times.txt").open("w") as f:
        for t in settle_times:
            f.write(f"{t / MYR_PER_CODE:.17g} 2\n")

    branch_times = [1, 5, 10, 20, 30, 40, 50, 60]
    with (HERE / "outputlists" / "branch_times.txt").open("w") as f:
        for t in branch_times:
            f.write(f"{t / MYR_PER_CODE:.17g} 2\n")

    settle = param
    settle = set_param(settle, "OutputDir", "outputs/settle")
    settle = set_param(settle, "OutputListFilename", "outputlists/settle_times.txt")
    settle = set_param(settle, "TimeMax", f"{SETTLE_MYR / MYR_PER_CODE:.17g}")
    settle = set_param(settle, "BHBenchmarkAccretionModel", "5")
    settle = set_param(settle, "BHBenchmarkStartTimeMyr", "1.0e9")
    (HERE / "param_settle.txt").write_text(settle)

    # Stage 2: branch the same settled snapshot into all six accretion laws.
    model_names = {
        0: "tng_bondi",
        1: "boosted_bondi",
        2: "am_bondi",
        3: "ffr_volume",
        4: "ffr_shell",
        5: "ffr_env",
    }
    params_dir = HERE / "params"
    params_dir.mkdir(parents=True, exist_ok=True)
    for model, name in model_names.items():
        branch = param
        branch = set_param(branch, "InitCondFile", "ics/settled_galaxy")
        branch = set_param(branch, "OutputDir", f"outputs/model_{model}_{name}")
        branch = set_param(branch, "OutputListFilename", "outputlists/branch_times.txt")
        branch = set_param(branch, "TimeMax", f"{(END_MYR - SETTLE_MYR) / MYR_PER_CODE:.17g}")
        branch = set_param(branch, "BHBenchmarkAccretionModel", str(model))
        branch = set_param(branch, "BHBenchmarkStartTimeMyr", "0.0")
        (params_dir / f"model_{model}_{name}.txt").write_text(branch)

    meta = {
        "seed": SEED,
        "box_kpc": BOX_KPC,
        "bh_mass_msun": M_BH,
        "settle_myr": SETTLE_MYR,
        "end_myr": END_MYR,
        "bh_accretion_radius_kpc": 0.30,
        "bh_feedback_radius_kpc": 0.60,
        "particle_counts": {"gas": n0, "dm": n1, "stars": n4, "bh": 1},
        "masses_msun": {
            "halo": M_HALO,
            "stellar_disk": M_STELLAR_DISK,
            "bulge": M_BULGE,
            "gas_extended": M_GAS_EXT,
            "gas_nuclear": M_GAS_NUC,
            "bh": M_BH,
        },
        "environment_model": {
            "thermal_floor": 0.4,
            "rotation_beta": 1.0,
            "description": "FFR_shell * Fth * Fwind * Frot; all support variables bounded in [0,1]",
        },
    }
    (HERE / "galaxy_test_metadata.json").write_text(json.dumps(meta, indent=2) + "\n")

    print(f"Wrote {path}")
    print(f"Particles: gas={n0}, DM={n1}, stars={n4}, BH=1; total={n0+n1+n4+1}")
    print(f"BH capture off until {SETTLE_MYR:g} Myr; run ends at {END_MYR:g} Myr")
    print("Wrote common settle run plus six post-settle accretion branches")


if __name__ == "__main__":
    main()
