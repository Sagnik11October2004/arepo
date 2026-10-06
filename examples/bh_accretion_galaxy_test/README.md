# Isolated galaxy BH accretion comparison

This test replaces repeated tuning against idealized one-flow boxes with a common evolving galactic environment.

## Setup

`prepare_galaxy_test.py` creates a deterministic HDF5 IC with a live dark halo, stellar disk plus bulge, extended gas disk plus a higher-resolution nuclear gas disk, and one central Type-5 black hole.

The BH contributes to gravity from the beginning. Resolved gas capture is disabled until 20 Myr via `BHBenchmarkStartTimeMyr`, allowing the galaxy to relax first. From 20 to 80 Myr, benchmark model 5 feeds the unresolved FFR-MACER reservoir. Feedback is disabled so the accretion estimators can be compared without a feedback loop changing the gas.

## Mode 5

The underlying law stays shell free-fall:

    Mdot_env = F_env Mdot_shell
    F_env = F_th F_wind F_rot

with bounded shell variables

    X_th   = c_s^2 / (c_s^2 + v_dyn^2)
    X_wind = v_coh^2 / (v_coh^2 + v_dyn^2 + c_s^2)
    X_rot  = v_phi^2 / (v_phi^2 + v_dyn^2 + c_s^2)

and

    F_th   = 1 - (1-f_th,min) X_th
    F_wind = (1+X_wind)^(-3/2)
    F_rot  = (1+beta_J X_rot^2)^(-1/2)

The exploratory defaults are `f_th,min=0.4` and `beta_J=1`. The correction is deliberately bounded and mild.

`v_coh` is the shell mass-weighted mean relative velocity multiplied by a velocity-coherence measure, so random turbulence is not treated as a coherent BHL wind.

## Run

    python3 prepare_galaxy_test.py
    ../bh_accretion_flow_tests/build_flow_tests.sh
    NTASKS=8 ./run_galaxy_test.sh

## Outputs

Every `BH_BENCHMARK_ALL` line contains the six raw estimators on the same evolving gas state: TNG Bondi, boosted Bondi-Hoyle, angular-momentum limited Bondi, volume FFR, shell FFR, and environment-corrected shell FFR.

The log also records `FENV`, `FTH`, `FWIND`, `FROT`, the bounded thermal/wind/rotation variables, velocity coherence, coherent Mach number, shell radius, shell sound speed, shell bulk speed, shell rotational speed, and shell dynamical speed.

`analyze_galaxy_accretion.py` writes:

- `galaxy_accretion_estimators.csv`
- `galaxy_bh_state.csv`
- `galaxy_accretion_rates.png`
- `galaxy_environment_factors.png`
- `galaxy_bh_reservoir_rates.png`

The BH-state CSV exports every available `BH_*` snapshot field, so reservoir mass, supply/feed/horizon/wind rates, accretion state, axes, energetic buffers, and benchmark realized rates remain available for later analysis.
