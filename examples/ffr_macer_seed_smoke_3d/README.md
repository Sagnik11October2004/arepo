# FFR-MACER FoF seeding smoke test

This is a runtime/conservation test for FFR-MACER FoF seeding, Iteration-6 inner accretion/energetics, Iteration-7 resolved bipolar wind feedback, Iteration-8 persistent-axis bipolar jet feedback, Iteration-9 feedback-aware BH timestep control, and Iteration-10 cached local DM velocity dispersion. It is not a production cosmology calculation.

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
2. every seed has a `BH_FFR: seed axis ... axis_source=local-gas|deterministic-fallback coherence=...` diagnostic consistent with `BHMinCoherence`; the local-gas estimate samples all resolved gas inside the proper `BHAccretionRadius` around the selected central cell before conversion;
3. Type-5 BHs in an HDF5 snapshot;
4. BH_Mass = 1e5 Msun;
5. BH_DiskMass = BH_MdotSupply = 0 because resolved capture is disabled;
6. N_gas + N_BH = 128^3 in every snapshot;
7. N_DM = 128^3 in every snapshot;
8. total snapshot mass conserved across the gas-to-BH conversion;
9. a dedicated checkpoint near a=0.02050 (z~47.7805) that contains no Type-5 BH and has z greater than the first seed redshift;
10. a snapshot after the first seeded snapshot in which the seeded BH ID still exists;
11. for the deliberately empty reservoir, BH_MdotFeed = BH_ProcessedEddRatio = 0, BH_MdotEdd is positive, and BH_Mode = 0 (ADIOS);
12. with an empty reservoir, BH_MdotHorizon, BH_MdotWind, BH_Lbol, BH_PWind, and BH_PJet remain exactly zero.

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

The runtime kernel checks each fired packet for exact returned wind mass, zero net bipolar kick momentum, and the requested kinetic-energy increment. Its startup self-test also checks that returned wind mass carrying the lobe-weighted ambient internal energy preserves the intended specific internal energy, including AREPO's cosmological `a^2` conserved-energy factor, before any mechanical kick. The verifier additionally requires both lobes to satisfy BHMinTargetsPerLobe. The helper raises BHJetBurstFactor to 1e30 so this remains an isolated wind regression after Iteration 8.

A copied checkpoint can be supplied with CHECKPOINT in the same way as the inner-flow helper.

## 9. Iteration-8 jet-axis and narrow-jet regression

The jet-only test again starts from the BH-free checkpoint. It suppresses wind release, uses Rfb=1.0 for adequate narrow-cone sampling in the coarse box, lowers the active target mass fraction to 0.01, and uses BHMinCoherence=0.10 so the persistent JetDir update is exercised:

    NTASKS=16 ./run_jet_restart.sh
    python3 verify_jet.py
    python3 inspect_bh.py output_jet

The verifier reads the per-step jet-axis logs and checks theta_new = theta_old exp(-dt/tdir), verifies unit DiscDir/JetDir vectors in snapshots, requires at least eight active targets in each narrow lobe, and checks every fired jet packet for exact kinetic-energy injection and zero net bipolar kick momentum. Jets return no separate rest mass; wind mass and energy are left buffered in this isolated test.

A copied checkpoint can be supplied with CHECKPOINT as for the other restart helpers.

## 10. Iteration-9 feedback/backlog timestep regression

The production limiter combines the nearest-gas hydro bin with the internal accuracy bound

    t_d = BHDiskTimeMyr * (1 + M_BH/M_d)^BHDiskTimeExponent
    dt_BH <= f_int min(t_d, Eth_w/P_w, Eth_j/P_j)

using physical time before AREPO performs its normal power-of-two timebin quantization. For a BH that completed an FFR transaction at the current synchronization point, the limiter reuses the exact beginning-of-transaction `t_d` already used by reservoir drainage and jet-axis relaxation; the timestep log records this as `tdFrozen=1`. Startup/new-seed/post-merger fallback states may recompute from the current state. If a surviving wind or jet buffer still contains at least BHMaxPacketsPerStep thresholds, it additionally forces at least one finer gravity bin, so a packet-cap backlog cannot wait on a coarse collisionless timestep.

The dedicated stress helper intentionally uses BHInternalTimestepFactor=4 and BHMaxPacketsPerStep=1. Those values are test-only: they make the ordinary accuracy limit loose enough that the independent backlog branch must activate. The backlog step is anchored to the gas+accuracy-limited candidate and may reduce it by only one additional factor of two per assignment; it is not recursively referenced to the already-shortened current BH bin. The same regression also checks that BH_MdotSupply is assembled as sum(dm_i/dt_i,hydro), avoiding diagnostic spikes when the BH timestep is shorter than the synchronized gas step. The run is shortened to a=0.02055 and still starts from the saved BH-free checkpoint:

    NTASKS=16 ./run_timestep_restart.sh
    python3 verify_timestep.py

Earlier Iteration-6/7/8 helpers set BHInternalTimestepFactor to a huge value internally so their historical regression sampling stays unchanged after Iteration 9.

## 11. Iteration-10 full-tree DM-dispersion regression

The production estimator runs only at valid full-gravity-tree synchronization points. By default it selects the globally nearest BHDMNeighbours Type-1 particles, computes the one-dimensional physical peculiar velocity dispersion about their mass-weighted local mean, and caches it in BH_SigmaDM between full-tree updates. Multimass zooms can override the compile-time BH_FFR_DM_TYPEMASK without changing All/BHP restart layouts. Transient imported gravity points carry velocity only when BLACKHOLE_FFR is compiled.

The isolated regression leaves BHFreeFallA=0 and disables the central binding term. Therefore the BH generates no wind/jet energy, while any positive burst threshold must come exclusively from SigmaDM:

    NTASKS=16 ./run_dm_restart.sh
    python3 verify_dm.py
    python3 inspect_bh.py output_dm

The verifier requires an exact 64-particle DM sample (no gas fallback), positive finite SigmaDM, and checks

    vbind^2 = SigmaDM^2
    EthWind = BHWindBurstFactor * 0.5 * Menc * vbind^2
    EthJet  = BHJetBurstFactor  * 0.5 * Menc * vbind^2

from runtime diagnostics. If fewer than BHDMNeighbours dark-matter particles exist, production code falls back to the nearest gas velocity sample and logs source=gas-fallback; no persistent fallback flag is added because the BHP restart ABI is frozen.

## 12. Iteration-11 dynamical-friction and simple-merger regressions

Iteration 11 adds a deliberately simple orbital sub-grid layer without changing the persistent BHP layout.

The dynamical-friction term reuses the exact nearest-DM sample from Iteration 10. The local mean and one-dimensional dispersion are mass-weighted, so multimass zoom particles do not bias the estimate. The transient DM environment is cached by BH ID, remains usable between full-tree refreshes while the BH has moved less than 0.25 Rfb, and is replicated across MPI ranks so domain migration does not silently disable drag. The gravity timestep also resolves BHInternalTimestepFactor times the cached Chandrasekhar damping time; the 50 percent kick cap remains only an emergency safeguard.

    NTASKS=16 ./run_df_restart.sh
    python3 verify_df.py

The DF and merger helpers use a second executable, `ArepoFFRSnapshotTest`, built automatically by `build_snapshot_ic.sh` from `Config_snapshot_ic.sh`. It intentionally omits `GENERATE_GAS_IN_ICS`, because these synthetic fresh-IC tests are constructed from AREPO hydrodynamic snapshots that already contain gas. The original `ArepoSeedTest` retains `GENERATE_GAS_IN_ICS` for the DM-only MUSIC initial conditions. The two executables use separate build directories so their compile-time configurations cannot contaminate one another.

For speed, the synthetic DF and merger helpers disable FoF seeding because their input files already contain the required BH(s). Their default evolution interval is `SHORT_DA=4e-6` (half the previous interval), and the test-only maximum-step fraction is `SHORT_STEP_FRAC=0.45` rather than `0.10`. This retains multiple synchronization points while avoiding unnecessary full mesh/gravity cycles. Both values can be overridden for debugging, for example `SHORT_DA=8e-6 SHORT_STEP_FRAC=0.1 ./run_merger_restart.sh`.

The merger rule is intentionally simple: on full synchronization points all live BHs are considered, while intermediate synchronization points additionally test the subset of BHs that are already gravity-active and drifted to the current time. Disjoint nearest pairs merge whenever their two proper accretion apertures overlap, d < 2 R_acc. The lowest particle ID survives, but the surviving particle is moved to the periodic dynamical-mass center of mass. Dynamical mass and linear momentum are conserved; BH mass, reservoir mass, wind mass/momentum/energy buffer, jet energy buffer and coherence are combined. Disc/jet directions fall back to the more massive progenitor when the merged directional vectors nearly cancel. Coalescence occurs after each progenitor's elapsed reservoir/DF update but before resolved wind/jet packets. The loser uses a temporary ID=0, Mass=0 Type-3 tombstone, which snapshot I/O excludes and the next normal domain rearrangement physically compacts.

The synthetic merger helper duplicates the one-BH Iteration-10 snapshot at a 0.05 proper separation and launches the result as a fresh HDF5 IC at that scale factor. This is deliberate: evolved Type-5 HDF5 snapshots are not restart-complete for FFR state and RestartFlag=2 now refuses them rather than silently destroying reservoir/buffer/directional state:

    NTASKS=16 ./run_merger_restart.sh
    python3 verify_merger.py
    python3 inspect_bh.py output_merger

The version-1 merger deliberately has no boundness, relative-velocity, binary-hardening, gravitational-wave delay or recoil criterion. Native RestartFlag=1 remains the supported continuation path for evolved FFR-MACER simulations; RestartFlag=2 is reserved for BH-free/pre-seed snapshots or deliberately fresh synthetic IC tests.

## 13. Iteration-11.5 hardening

Iteration 11.5 is a correctness/robustness pass rather than a new accretion model. It rejects evolved Type-5 RestartFlag=2 snapshots, versions the native FFR restart payload, serializes the rank-local dynamical-friction cache separately from the frozen BHP record and re-replicates it after RestartFlag=1, bounds DF-cache reuse by both 0.25 Rfb displacement and a physical age limit, preserves merger timestep/backlog information, tightens code-unit conservation tolerances, and logs the unresolved kinetic-energy dissipation of inelastic gas capture.

Feedback overlap conflicts still permit at most one event to own any active gas cell in a packet round, but priority is now a deterministic function of BH ID, synchronization time and packet round instead of permanently favoring the lower ID. This removes systematic starvation while remaining MPI-order independent.

For multimass zoom ICs, set BH_FFR_DM_TYPEMASK at compile time to the sum of the selected particle-type bits. Type 1 remains the default. For example, Types 1, 2 and 3 use:

    BH_FFR_DM_TYPEMASK=2+4+8

The nearest-neighbour communication payload now defaults to 64 entries, which covers the fiducial 40-neighbour estimator and the 64-neighbour smoke test without carrying the previous 256-candidate result structure for every exported BH. Runs that deliberately need a larger sample can set `BH_FFR_DM_MAX_NEIGHBOURS=<N>` in Config.sh; parameter validation uses the compiled capacity.

The current native restart format is FFR restart version 2. It stores the frozen BHP record size and the separate DF-cache record size and fails explicitly on an incompatible layout rather than interpreting shifted bytes as valid state. Restart files written by older FFR binaries should be resumed with the binary that wrote them.

After pulling the hardening branch, rerun the regression chain before coupled validation. The regression runners automatically invoke build.sh if ArepoSeedTest is missing:

    NTASKS=16 ./run_dm_restart.sh && python3 verify_dm.py
    NTASKS=16 ./run_df_restart.sh && python3 verify_df.py
    NTASKS=16 ./run_merger_restart.sh && python3 verify_merger.py
    NTASKS=4  ./run_restart_guard.sh

The merger verifier accepts both explicit per-particle Masses datasets and the standard Gadget/AREPO Header/MassTable convention for fixed-mass particle species.

## Why the boosted sigma8?

The purpose here is to exercise MPI FoF seeding and conservative gas-to-BH conversion quickly in a tiny box. sigma8=2.0 is deliberately non-production. After this passes, repeat at 256^3 with the intended physical cosmology/threshold before treating the seeding prescription as scientifically validated.
