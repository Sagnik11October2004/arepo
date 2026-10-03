# FFR-MACER FoF seeding smoke test

This is a runtime/conservation test for FFR-MACER FoF seeding plus the Iteration-5 empty-reservoir state selector. It is not a production cosmology calculation.

## Test design

- periodic box: 1 cMpc/h
- initial redshift: z=49
- explicit pre-seed checkpoint: a=0.02050 (z~47.7805), before the deterministic first seed near a~0.020526 (z~47.72)
- post-seed snapshot: a=0.02056 (z~47.64)
- final scale factor: a=0.02060 (z~47.54), a short persistence interval after the first BH snapshot
- unigrid resolution: 128^3 parent DM particles
- AREPO GENERATE_GAS_IN_ICS splits these into 128^3 gas cells + 128^3 Type-1 DM particles
- Omega_m=0.31, Omega_b=0.048, Omega_Lambda=0.69, h=0.68, n_s=0.96
- sigma8=2.0 **only for this smoke test**, to promote early halo formation before z=20
- FoF seed threshold: 2e6 Msun
- BH seed mass: 1e5 Msun
- seed redshift condition: z>20
- maximum donor removal fraction: 0.5
- resolved post-seed gas capture is disabled with BHFreeFallA=0

The expected mass resolution is approximately 5.10e4 Msun per Type-1 DM particle and 9.34e3 Msun per gas cell. A 2e6 Msun halo is therefore close to the minimum resolved FoF scale and is intentionally chosen to make this a cheap software test.

## 1. Pull the test

    cd ~/Research/arepo
    git checkout ffr-macer
    git pull

## 2. Build MUSIC2

The current MUSIC2 repository supports AREPO output, 2LPT, FFTW3, GSL and HDF5. The helper defaults to ~/Research/MUSIC2:

    cd ~/Research/arepo/examples/ffr_macer_seed_smoke_3d
    chmod +x install_music2.sh generate_ics.sh build.sh run.sh
    ./install_music2.sh

If MUSIC2 is already installed, skip this step and set MUSIC_BIN when generating the IC.

## 3. Generate the IC

    ./generate_ics.sh

or, for a non-default MUSIC executable,

    MUSIC_BIN=/absolute/path/to/MUSIC ./generate_ics.sh

The generated file is ics.hdf5. MUSIC2 writes a DM-only AREPO HDF5 unigrid; AREPO then performs the baryon split at startup using GENERATE_GAS_IN_ICS.

## 4. Build the dedicated AREPO executable

    ./build.sh

This builds ../../ArepoSeedTest with this directory's Config.sh. It does not alter Config_FFR_MACER.sh.

## 5. Run with MPI

Default: four MPI ranks. The checked-in 1500 MB/rank memory cap also permits a 16-rank smoke run on a ~32 GB machine.

    ./run.sh

To change the rank count:

    NTASKS=16 ./run.sh

The run writes output/ and run.log.

## 6. Verify

    python3 verify_seed.py

A passing test requires all of the following:

1. at least one "BH_FFR: seeded FoF group" event at z>20;
2. Type-5 BHs in an HDF5 snapshot;
3. BH_Mass = 1e5 Msun;
4. BH_DiskMass = BH_MdotSupply = 0 because resolved capture is disabled;
5. N_gas + N_BH = 128^3 in every snapshot;
6. N_DM = 128^3 in every snapshot;
7. total snapshot mass conserved across the gas-to-BH conversion;
8. a dedicated checkpoint near a=0.02050 (z~47.7805) that contains no Type-5 BH and has z greater than the first seed redshift;
9. a snapshot after the first seeded snapshot in which the seeded BH ID still exists;
10. for the deliberately empty reservoir, BH_MdotFeed = BH_ProcessedEddRatio = 0, BH_MdotEdd is positive, and BH_Mode = 0 (ADIOS);
11. BH_MdotHorizon, BH_MdotWind, BH_Lbol, BH_PWind, and BH_PJet remain zero, confirming that Iteration-6+ physics has not leaked into Iteration 5.

If the deterministic test realization does not seed near z~47.72, that is a test failure or a changed realization/configuration. Do not lower the threshold silently; inspect run.log first.

## Why the boosted sigma8?

The purpose here is to exercise MPI FoF seeding and conservative gas-to-BH conversion quickly in a tiny box. sigma8=2.0 is deliberately non-production. After this passes, repeat at 256^3 with the intended physical cosmology/threshold before treating the seeding prescription as scientifically validated.
