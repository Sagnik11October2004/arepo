# Isolated galaxy BH accretion comparison

This is the main development test for the accretion prescriptions. It uses one common evolving galaxy rather than tuning mode 5 against the idealized 42-flow matrix.

## Experiment design

1. `prepare_galaxy_test.py` creates a deterministic live galaxy: Hernquist-like DM halo, stellar disk + bulge, extended gas disk + resolved nuclear gas disk, and one central Type-5 BH.
2. The BH contributes to gravity but gas capture is disabled for 20 Myr.
3. The final settled snapshot is converted into a fresh HDF5 IC with the clock and BH subgrid state reset.
4. Six branches start from that exact same settled phase-space state, one for each accretion model, and evolve for another 60 Myr with the same FFR-MACER reservoir/inner-flow machinery and feedback disabled.

This isolates the accretion prescription itself while giving every model the same relaxed galaxy as its starting point.

## Mode 5

The underlying law remains shell free-fall:

    Mdot_env = F_env Mdot_shell
    F_env = F_th F_wind F_rot

with bounded support variables

    X_th   = c_s^2 / (c_s^2 + v_dyn^2)
    X_wind = v_coh^2 / (v_coh^2 + v_dyn^2 + c_s^2)
    X_rot  = v_phi^2 / (v_phi^2 + v_dyn^2 + c_s^2)

and mild corrections

    F_th   = 1 - (1-f_th,min) X_th
    F_wind = (1+X_wind)^(-3/2)
    F_rot  = (1+beta_J X_rot^2)^(-1/2)

The current exploratory values are `f_th,min=0.4` and `beta_J=1`. Because all X variables are bounded in [0,1], this cannot turn mode 5 into an arbitrarily strong Bondi-like suppression of a low-mass seed.

`v_coh` is the shell mass-weighted mean relative velocity multiplied by a velocity-coherence measure, so random turbulence cancels rather than being treated as a coherent wind.

## Six post-settle branches

- 0: TNG Bondi
- 1: boosted Bondi-Hoyle
- 2: angular-momentum limited Bondi
- 3: volume FFR
- 4: shell FFR
- 5: environment-corrected shell FFR

All six feed the same unresolved reservoir (`BHBenchmarkAccretionTarget=0`) and run with feedback disabled. The `BH_BENCHMARK_ALL` line still evaluates all six raw estimators in every branch, so both counterfactual instantaneous rates and true divergent evolutionary histories are available.

## Run

    cd ~/Research/arepo
    git pull --ff-only origin ffr-macer-convergence-benchmarks
    cd examples/bh_accretion_flow_tests
    ./build_flow_tests.sh
    cd ../bh_accretion_galaxy_test
    NTASKS=8 ./run_galaxy_test.sh

`run_galaxy_test.sh` automatically prepares the initial galaxy, runs the 20 Myr settling stage, creates the fresh settled branch IC, runs all six 60 Myr branches, and launches the analysis.

## Diagnostics

`BH_BENCHMARK_ALL` records time, whether capture is enabled, all six raw estimators, and the complete mode-5 environmental state: `FENV`, `FTH`, `FWIND`, `FROT`, bounded thermal/wind/rotation variables, velocity coherence, coherent Mach number, shell radius, sound speed, bulk speed, rotational speed, and dynamical speed.

`analyze_galaxy_accretion.py` writes:

- `galaxy_accretion_estimators.csv` — all raw estimator samples in the settle and six branch runs;
- `galaxy_bh_state.csv` — every available `BH_*` snapshot field, flattened into a table;
- `galaxy_branch_selected_rates.png` — selected supply law in each evolving branch;
- `galaxy_branch_bh_mass.png` — BH mass histories;
- `galaxy_environment_factors.png` — the mode-5 correction factors.

The CSV state table includes reservoir mass, supply/feed/horizon/wind rates, Eddington ratio, accretion mode, orientation axes, feedback buffers, benchmark realized rate, and any other `BH_*` field present in the snapshot.
