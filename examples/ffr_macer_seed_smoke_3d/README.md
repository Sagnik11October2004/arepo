# FFR-MACER FoF seeding smoke test

This is a runtime/conservation test for FFR-MACER FoF seeding, Iteration-6 inner accretion/energetics, Iteration-7 resolved bipolar wind feedback, and Iteration-8 persistent-axis bipolar jet feedback. It is not a production cosmology calculation.

## Test design

- periodic box: 1 cMpc/h
- initial redshift: z=49
- explicit pre-seed checkpoint: a=0.02050 (z~47.7805), before the deterministic first seed near a~0.020526 (z~47.72)
- post-seed snapshot: a=0.02056 (z~47.64)
- final scale factor: a=0.02060 (z~47.54), a short persistence interval after the first BH snapshot
- maximum timestep: 5e-4 in scale factor/log-a control, below the shortened total interval TimeMax-TimeBegin=6e-4 and much less restrictive than the previous 1e-4 smoke-test cap
- unigrid resolution: 128^3 parent DM particles
- AREPO GENERATE_GAS_IN_ICS splits these into 128^3 gas cells + 128^3 Type-1 DM particles
- Omega_m=0.31, Omega_b=0.048, Omega_Lambda=0.69, h=0.68, n_s=0.96
- sigma8=2.0 **only for this smoke test**, to promote early halo formation before z=20
- FoF seed threshold: 2e6 Msun
- BH seed mass: 1e5 Msun
- seed redshift condition: z>20
- maximum donor removal fraction: 0.5
- the canonical zero-capture regression keeps BHFreeFallA=0
- the Iteration-6 inner-flow regression restarts from the BH-free checkpoint and enables capture with BHFreeFallA=1e-3

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
11. with an empty reservoir, BH_MdotHorizon, BH_MdotWind, BH_Lbol, BH_PWind, and BH_PJet remain exactly zero.

If the deterministic test realization does not seed near z~47.72, that is a test failure or a changed realization/configuration. Do not lower the threshold silently; inspect run.log first.

## 7. Iteration-6 inner-flow regression

After the canonical run has produced the BH-free `output/snap_001.hdf5`, rebuild with the current branch and run:

    NTASKS=16 ./run_inner_restart.sh
    python3 verify_inner.py
    python3 inspect_bh.py output_inner

The helper copies the BH-free checkpoint into `output_inner/`, sets `BHFreeFallA=1e-3`, and restarts with RestartFlag=2 from snapshot 001. It therefore begins near z=47.79 rather than repeating the z=49 evolution.

The Iteration-6 verifier checks the conservative inner-flow identity

    BH_MdotFeed = BH_MdotHorizon + BH_MdotWind

along with the dynamical-mass ledger, hot-flow retention, hot wind/jet powers, radiative luminosity, cold Gofford closure when sampled, and monotonic accumulation of wind/jet buffers. The helper raises the wind/jet burst factors to 1e30 so this regression remains an isolated Iteration-6 test even after feedback is implemented.

If the checkpoint was copied elsewhere, use for example:

    CHECKPOINT=checkpoints/preseed_a0.0204956088_z47.790939.hdf5 NTASKS=16 ./run_inner_restart.sh

## 8. Iteration-7 bipolar wind-feedback regression

The wind-feedback test uses the same BH-free checkpoint, enables capture, enlarges the feedback aperture to 1 proper code-length unit so both lobes contain many cells in this coarse smoke box, and lowers the wind burst factor to exercise actual packet release:

    NTASKS=16 ./run_feedback_restart.sh
    python3 verify_feedback.py
    python3 inspect_bh.py output_feedback

The runtime kernel checks each fired packet for exact returned wind mass, zero net bipolar kick momentum, and the requested kinetic-energy increment. The verifier additionally requires both lobes to satisfy BHMinTargetsPerLobe. The helper raises BHJetBurstFactor to 1e30 so this remains an isolated wind regression after Iteration 8.

A copied checkpoint can be supplied with CHECKPOINT in the same way as the inner-flow helper.

## 9. Iteration-8 jet-axis and narrow-jet regression

The jet-only test again starts from the BH-free checkpoint. It suppresses wind release, uses Rfb=1.0 for adequate narrow-cone sampling in the coarse box, lowers the active target mass fraction to 0.01, and uses BHMinCoherence=0.10 so the persistent JetDir update is exercised:

    NTASKS=16 ./run_jet_restart.sh
    python3 verify_jet.py
    python3 inspect_bh.py output_jet

The verifier reads the per-step jet-axis logs and checks theta_new = theta_old exp(-dt/tdir), verifies unit DiscDir/JetDir vectors in snapshots, requires at least eight active targets in each narrow lobe, and checks every fired jet packet for exact kinetic-energy injection and zero net bipolar kick momentum. Jets return no separate rest mass; wind mass and energy are left buffered in this isolated test.

A copied checkpoint can be supplied with CHECKPOINT as for the other restart helpers.

## Why the boosted sigma8?

The purpose here is to exercise MPI FoF seeding and conservative gas-to-BH conversion quickly in a tiny box. sigma8=2.0 is deliberately non-production. After this passes, repeat at 256^3 with the intended physical cosmology/threshold before treating the seeding prescription as scientifically validated.
