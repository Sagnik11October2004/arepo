# High-redshift rare-peak BH growth experiment

This example builds a controlled cosmological rare-peak experiment for comparing
black-hole accretion and feedback prescriptions on the same forming halo.

It is **not** a representative cosmological volume: the box deliberately uses
`sigma_8=1.2` to modestly enhance early structure in a small volume.

## Design

- box: 0.5 cMpc/h
- cosmology: Omega_m=0.31, Omega_b=0.048, Omega_Lambda=0.69, h=0.68, n_s=0.96
- moderately boosted amplitude: sigma8=1.2
- MUSIC start: z=49, 2LPT
- quick-scout parent unigrid: 64^3 DM particles
- AREPO startup split: 64^3 gas + 64^3 Type-1 DM
- expected gas mass: about 9.34e3 Msun
- expected DM mass: about 5.10e4 Msun
- simple primordial AREPO cooling enabled; star formation disabled
- FoF compiled throughout
- seed mass: 1e5 Msun
- eligible host threshold: 5e6 Msun
- dedicated compile-time single-seed mode chooses the globally most massive
  eligible FoF halo and prevents later creation of additional BH seeds
- seed window: 20 < z < 22
- gas capture and BH feedback are disabled until the common z=20 branch state

The seed window is enforced in the source, not by a runner-time hack.  Both
rare-peak executables define `BH_FFR_SEED_MAX_REDSHIFT=22`, while the runtime
parameter remains `BHSeedMinRedshift=20`.  The FoF seed transaction therefore
requires the strict condition

```text
20 < z < 22
```

before any eligible halo can be seeded.  The upper edge is compile-time gated
rather than added to `All`, so the existing BLACKHOLE_FFR restart layout and
the default one-sided production seeding behaviour are unchanged.

The run is still split into two stages for reproducibility:

1. z=49 -> z=22, which must remain BH-free because the strict upper edge is not
   yet crossed.
2. z=22 -> z=20, where the single-most-massive eligible FoF halo may be seeded.

The second stage starts from the exact BH-free z=22 hydrodynamic snapshot using
the executable that omits `GENERATE_GAS_IN_ICS`.

## Files

- `music.conf` -- z=49, 64^3 MUSIC2 scout IC configuration.
- `Config_initial.sh` -- initial DM-only MUSIC input; includes gas generation.
- `Config_evolved.sh` -- already-split gas+DM cosmology and science branches.
- `param_preseed.txt` -- z=49 -> z=22, no seeding.
- `param_seed_window.txt` -- z=22 -> z=20, one seed, no capture/feedback.
- `run_to_z20.sh` -- complete common-state runner.
- `prepare_branch_params.py` -- writes all accretion/feedback branch parameters.
- `run_branches.sh` -- runs a selectable branch matrix.

## 1. Generate the IC

The helper assumes MUSIC2 is at `~/Research/MUSIC2/build/MUSIC`.

```bash
cd ~/Research/arepo
git checkout ffr-macer-convergence-benchmarks
git pull --ff-only origin ffr-macer-convergence-benchmarks

cd examples/bh_cosmo_rare_peak_3d
chmod +x generate_ics.sh build.sh run_to_z20.sh run_branches.sh

./generate_ics.sh
```

For another MUSIC executable:

```bash
MUSIC_BIN=/absolute/path/to/MUSIC ./generate_ics.sh
```

`verify_ic.py` checks z=49, the 64^3 parent count, and that the expected
post-split gas mass is below 1e4 Msun. `generate_ics.sh` and `run_to_z20.sh`
also print the actual free space on the filesystem and refuse to proceed below
5 GiB by default (`MIN_FREE_GIB` can override the guard).

## 2. Build

```bash
./build.sh
```

This produces two independent executables/build directories:

- `ArepoRarePeakInitial`: includes `GENERATE_GAS_IN_ICS`.
- `ArepoRarePeak`: reads already-split cosmological snapshots/branch ICs.

Both configurations include `INPUT_IN_DOUBLEPRECISION`,
`OUTPUT_IN_DOUBLEPRECISION`, HDF5, cooling, FoF and BLACKHOLE_FFR.

## 3. Evolve to the common z=20 seeded state

```bash
NTASKS=4 ./run_to_z20.sh
```

If Stage A has already produced one or more `output_preseed/snap_*.hdf5`
snapshots, it can be resumed from the latest existing snapshot after rebuilding
the evolved executable:

```bash
RESUME_PRESEED=1 NTASKS=4 ./run_to_z20.sh
```

The resume path validates that the selected checkpoint is BH-free, uses AREPO
`RestartFlag=2`, and runs it with `ArepoRarePeak` (the executable without
`GENERATE_GAS_IN_ICS`). It does not delete the existing Stage-A snapshots.


The rare-peak executables also enable `HIERARCHICAL_GRAVITY` and `ALLOW_DIRECT_SUMMATION` (`DIRECT_SUMMATION_THRESHOLD=500`). This matters once only a few particles occupy the shortest time bins: the non-hierarchical gravity path otherwise rebuilds a full multi-million-particle tree at every tiny synchronization step. `CELL_CENTER_GRAVITY` follows the public AREPO cosmological example. These are gravity-integration/performance choices; they do not change the BH accretion prescription.

The checked-in rare-peak parameter files use `MaxMemSize=1500` MB per MPI
rank.  This is intentionally a per-rank cap: for example, 4 ranks advertise
at most 6 GB to AREPO's startup memory guard, appropriate for the 64^3 scout.  Do not interpret `MaxMemSize` as a node-wide total.

The runner:

1. evolves z=49 -> z=22 with no BH;
2. verifies the last stage-A snapshot is BH-free and at z=22;
3. copies that state to `ics/z22_preseed.hdf5`;
4. evolves z=22 -> z=20 with FoF seeding active;
5. audits the FoF seed log and requires exactly one event with `20<z<22`;
6. records seed redshift, FoF host mass, FoF group and BH ID in
   `ics/common_seed_metadata.json`;
7. requires exactly one 1e5 Msun BH at z=20 with the same BH ID;
8. writes `ics/common_z20_seeded.hdf5`.

Because capture is disabled through this stage, the common BH must still have
exactly the seed mass. If no eligible >=5e6 Msun halo appears in the requested
window, the runner exits rather than silently lowering the host threshold.

## 4. Accretion modes available after z=20

All six current estimators are exposed:

| ID | Model |
|---:|---|
| 0 | TNG Bondi |
| 1 | boosted Bondi-Hoyle |
| 2 | angular-momentum-limited Bondi |
| 3 | volume FFR |
| 4 | shell FFR |
| 5 | convergence/J-corrected shell FFR |

All three feedback choices are also generated:

- `none`: feedback off; direct target by default, selectable with
  `--none-target reservoir`.
- `tng`: TNG thermal/kinetic feedback, forced to the direct target.
- `macer`: MACER wind/jet feedback, forced to the reservoir target.

Generate the full parameter matrix without running it:

```bash
python3 prepare_branch_params.py --z-end 15
```

This writes 18 files under `params/`.

## 5. Run selected branches

Do not launch all 18 on a laptop first. A useful first cosmological check is:

```bash
MODELS="0 4 5" FEEDBACKS="none" ZEND=19 NTASKS=8 ./run_branches.sh
```

A later full matrix is:

```bash
MODELS="0 1 2 3 4 5" FEEDBACKS="none tng macer" ZEND=15 NTASKS=8 ./run_branches.sh
```

Outputs are kept separate as
`outputs/<accretion-model>__<feedback>/`.

## Current mode-5 normalization

Mode 5 is now interpreted as a resolved dynamical suppression of the
geometric shell free-fall supply,

```text
Mdot_5 = Mdot_shell,geom * f_M * f_j
f_j = [1 + (r_circ/R_acc)^2]^(-1/2)
```

with no additional `1e-3` normalization.  Here
`r_circ/R_acc = j^2/(G M_cen R_acc) = v_phi^2/v_K^2` is also the ratio of
centrifugal to gravitational acceleration at the aperture.  The adopted
`f_j` is a smooth closure: it approaches unity for negligible rotation and
the gravity/centrifugal-force ratio for strong rotational support.  The
specific quadrature form is a model choice rather than a unique first-principles
derivation and should be tested against resolved inflow calculations.

If a future residual global efficiency is introduced for mode 5, it should
represent missing unresolved physics and be restricted to the range 0.1--1;
it must not reintroduce the retired `1e-3` prefactor.


## Native-restart z~10 science suite

The current controlled comparison starts from the untouched native AREPO
restart written immediately after the first seed event:

- z_seed = 20.8837
- BH ID = 1000120011
- M_BH = 1e5 Msun
- four MPI ranks in the checkpoint
- no resolved gas capture had occurred before the checkpoint

BH-containing FFR/MACER states must be branched with native
`RestartFlag=1`; a snapshot restart is not state-safe because the unresolved
reservoir, feedback buffers, directions, and timestep bookkeeping live in the
serialized BH state.

The checked-in branch parameter files are under `params_native/`:

- `ffr_shell_A1e2__macer.txt`: shell FFR with A_ff=1e-2, reservoir, MACER.
- `ffr_shell_A1e1__macer.txt`: shell FFR with A_ff=1e-1, reservoir, MACER.
- `convj_shell_ffr__macer.txt`: geometric shell x f_M x f_j, no extra
  suppression, reservoir, MACER.
- `tng_bondi__tng.txt`: TNG Bondi, direct target, TNG feedback.
- `ffr_shell_A1e2__none.txt`, `convj_shell_ffr__none.txt`, and
  `tng_bondi__none.txt`: feedback-free controls.

They share `output_branch_z10.txt` and request `TimeMax=0.09`, which is
chosen for the native integer timeline and corresponds to an actual endpoint
near z=10.11 for this checkpoint.  The branch files use the standard AREPO
`MinSizeTimestep=1e-7`.  No additional BLACKHOLE_FFR minimum-timestep floor
is added; if the resolved hydro calculation genuinely requires a smaller
standard AREPO timestep, that should be diagnosed rather than hidden by a
BH-specific clamp.

The runner is deliberately non-destructive and requires an explicit branch
selection:

```bash
chmod +x run_native_restart_suite.sh

BRANCHES="convj_shell_ffr__macer" NTASKS=4 ./run_native_restart_suite.sh
```

For several branches, list them explicitly:

```bash
BRANCHES="ffr_shell_A1e2__macer ffr_shell_A1e1__macer convj_shell_ffr__macer tng_bondi__tng" \
  NTASKS=4 ./run_native_restart_suite.sh
```

By default it reads the pristine checkpoint from
`restart_archive/common_seed_spare/restartfiles`, copies it independently
for each selected branch, and refuses to overwrite an existing output or log.
The old `run_branches.sh` snapshot workflow is disabled.

## Scientific interpretation

This box should be described as a controlled overdense/rare-peak realization.
It is intended to answer a differential question:

> Starting from the same halo, gas distribution and 1e5 Msun seed, how do
> different unresolved accretion and feedback prescriptions change early BH
> growth and the nuclear gas supply?

It should not be used to infer the cosmological abundance of such systems.

The checked-in cooling model is AREPO's simple primordial cooling, not the
non-equilibrium H2/HD chemistry used in the separate primordial-collapse work.
