# High-redshift rare-peak BH growth experiment

This example builds a controlled cosmological rare-peak experiment for comparing
black-hole accretion and feedback prescriptions on the same forming halo.

It is **not** a representative cosmological volume: the box deliberately uses
`sigma_8=2` to produce an early massive halo in a small volume.

## Design

- box: 0.5 cMpc/h
- cosmology: Omega_m=0.31, Omega_b=0.048, Omega_Lambda=0.69, h=0.68, n_s=0.96
- deliberately boosted amplitude: sigma8=2
- MUSIC start: z=99, 2LPT
- parent unigrid: 128^3 DM particles
- AREPO startup split: 128^3 gas + 128^3 Type-1 DM
- expected gas mass: about 1.17e3 Msun
- expected DM mass: about 6.37e3 Msun
- simple primordial AREPO cooling enabled; star formation disabled
- FoF compiled throughout
- seed mass: 1e5 Msun
- eligible host threshold: 1e7 Msun
- dedicated compile-time single-seed mode chooses the globally most massive
  eligible FoF halo and prevents later creation of additional BH seeds
- seed window: 20 < z < 22
- gas capture and BH feedback are disabled until the common z=20 branch state

The seed window is implemented as two stages rather than changing the validated
runtime redshift gate:

1. z=99 -> z=22 with `BHSeedMinRedshift=100`, so seeding is disabled.
2. z=22 -> z=20 with `BHSeedMinRedshift=20`.

The second stage starts from the exact BH-free z=22 hydrodynamic snapshot using
the executable that omits `GENERATE_GAS_IN_ICS`.

## Files

- `music.conf` -- z=99 MUSIC2 IC configuration.
- `Config_initial.sh` -- initial DM-only MUSIC input; includes gas generation.
- `Config_evolved.sh` -- already-split gas+DM cosmology and science branches.
- `param_preseed.txt` -- z=99 -> z=22, no seeding.
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

`verify_ic.py` checks z=99, the 128^3 parent count, and that the expected
post-split gas mass is below 1e4 Msun.

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
NTASKS=8 ./run_to_z20.sh
```

The runner:

1. evolves z=99 -> z=22 with no BH;
2. verifies the last stage-A snapshot is BH-free and at z=22;
3. copies that state to `ics/z22_preseed.hdf5`;
4. evolves z=22 -> z=20 with FoF seeding active;
5. requires exactly one 1e5 Msun BH;
6. writes `ics/common_z20_seeded.hdf5`.

Because capture is disabled through this stage, the common BH must still have
exactly the seed mass. If no eligible >=1e7 Msun halo appears in the requested
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

## Important mode-5 normalization note

The current development implementation of mode 5 deliberately uses the
unit-efficiency geometric shell rate times `FMACH*FJ` so its *shape* can be
tested cleanly against the idealized flow suite. It does not yet apply the
production `A_ff=1e-3` normalization used by the original shell FFR.

Therefore mode 5 is wired into the cosmological branch runner, but its final
`A_ff`/free-fall normalization must be frozen before the mode-5 branch is
interpreted as a production science result.

## Scientific interpretation

This box should be described as a controlled overdense/rare-peak realization.
It is intended to answer a differential question:

> Starting from the same halo, gas distribution and 1e5 Msun seed, how do
> different unresolved accretion and feedback prescriptions change early BH
> growth and the nuclear gas supply?

It should not be used to infer the cosmological abundance of such systems.

The checked-in cooling model is AREPO's simple primordial cooling, not the
non-equilibrium H2/HD chemistry used in the separate primordial-collapse work.
