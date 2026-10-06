# Black-hole accretion flow IC suite

This directory generates the controlled non-spherical IC hierarchy used after the completed spherical free-fall benchmark.

The generator creates exactly 42 double-precision HDF5 ICs:

- BHL/moving flow: Mach 0, 0.5, 1, 2, 5 at L0/L1/L2 (15 files).
- Rotating inflow: f_rot = 0, 0.25, 0.5, 0.75, 1 at L0/L1/L2 (15 files).
- Turbulent inflow: turbulent Mach 0.5, 1, 2, 5 at L0/L1/L2 (12 files).

Common geometry: 100 pc box, central 1e5 Msun Type-5 BH, R_acc=12.5 pc, R_feedback=25 pc. L0/L1/L2 are nested 16^3/32^3/64^3 gas hierarchies. Each split conserves each parent mass and vector momentum to floating-point precision.

Run:

```bash
python3 generate_all_ics.py
python3 verify_ics.py
```

The generated HDF5 files and per-IC metadata are intentionally git-ignored. `ic_manifest.json` and `verification_report.json` are also regenerated locally.

## BHL

The upstream adiabatic sound speed is 5 km/s. The uniform density is computed, not hard-coded, so that the analytic BHL rate at Mach 1 is 1e-4 Msun/yr. Gas flows in +x past a stationary BH. Metadata records the analytic BHL rate and focusing radius for every Mach number and resolution.

## Rotating inflow

This starts from the spherical inflow construction and adds coherent rotation about z:

v_phi = f_rot sqrt(G M_BH / r_eff)

v_r = -sqrt((2-f_rot^2) G M_BH / r_eff).

The density is recomputed from the fixed radial supply 1e-4 Msun/yr, so changing f_rot does not simply change the imposed outer mass supply. The inner profile is regularized at 0.5 pc; flow is unmodified through 35 pc and smoothly tapered to a quiet 1%-density ambient by 40 pc.

## Turbulent inflow

This uses the f_rot=0 spherical inflow plus a fixed-seed solenoidal Fourier field. Only integer modes with |k|<=4 are used, corresponding to wavelengths >=25 pc in the 100 pc box, so all modes are resolved already at L0. Vector amplitudes are projected perpendicular to k and scaled as |k|^-2, giving P_v(k) proportional to k^-4. A single common realization is evaluated at every resolution. Bulk momentum and coherent rigid rotation are removed from the continuous field using the L0 mass distribution; the amplitude is normalized inside 25 pc. Fine-level child velocities then receive only the common per-parent momentum correction required by the nested hierarchy.

The verifier checks HDF5 structure, particle counts, BH state, positivity/finiteness, shell populations, profile/Mach targets, and exact L0->L1->L2 parent-child mass and momentum closure.

## Build and run matrix

The canonical compile configuration is `Config_flow_tests.sh`. It reads and writes
double-precision HDF5 data and includes:

```text
DOUBLEPRECISION=1
DOUBLEPRECISION_FFTW
INPUT_IN_DOUBLEPRECISION
OUTPUT_IN_DOUBLEPRECISION
OUTPUT_COORDINATES_IN_DOUBLEPRECISION
NGB_TREE_DOUBLEPRECISION
```

Build the dedicated executable with:

```bash
./build_flow_tests.sh
```

Generate the parameter files and output lists with:

```bash
python3 generate_run_files.py
python3 verify_run_files.py
```

The generated matrix contains 60 planned passive runs:

- 42 one-kyr startup diagnostics: every BHL, rotating and turbulent IC at L0/L1/L2.
- 18 representative evolved-flow runs at L1/L2:
  - BHL Mach 0.5, 2 and 5, each evolved for one local `R_acc/v_eff` time.
  - rotating `f_rot=0.25,0.75,1`, each evolved for 200 kyr.
  - turbulent Mach 0.5, 2 and 5, each evolved for 200 kyr.

All of these use direct accretion, no feedback and
`BHBenchmarkEddingtonFactor=1e-8`, so the selected shell-FFR sink is dynamically
negligible. The diagnostic line `BH_BENCHMARK_ALL` reports TNG Bondi, boosted
BHL, AM-Bondi, volume FFR, shell FFR, and the cell-wise effective-infall FFR raw
rate from the same gas state. Benchmark model 5 keeps the shell/free-fall
backbone and changes only the cell fall time:

`Mdot_5 = C_shell * sum_i m_i f_i / t_ff,i`

with

`f_i = [1 + (c_s^2 + v_perp^2 + v_r,out^2)/(v_g^2 + v_r,in^2)]^(-3/2)`,
`v_g^2 = G M_central/r_i`,
`v_r,in=max(-v_r,0)`, and `v_r,out=max(v_r,0)`.

Thus thermal, transverse/rotational, and outward radial motion lengthen the
fall time, while already-inward radial motion reduces that delay.  There is
no aperture-level Mach switch, Bondi cap, or extra M_dyn term.  The diagnostic
line prints `FFR_SHELL_GEOM`, `FFR_EFFECTIVE`, and the shell-rate-weighted
mean `effsupport`.

Every output-list entry uses `DumpFlag=2`: a full snapshot is written but the
compiled FoF module is not run. FoF remains in the build only because the
production BLACKHOLE_FFR seeding implementation requires it at compile time.

Run a complete stage only after the single-case smoke test has passed:

```bash
./run_flow_stage.sh diagnostic
./run_flow_stage.sh evolved
```

The stage runner defaults to 2/8/16 MPI ranks for L0/L1/L2. Set
`EXEC`, `MPIEXEC`, or `KEEP_RESTARTS` in the environment to override the
executable, MPI launcher, or final-restart cleanup policy. Active,
model-specific sink/feedback evolutions are intentionally not generated yet;
their matrix should be chosen only after the passive diagnostic/evolved-flow
results are inspected.

### Case-appropriate reference comparison

Use:

```bash
python3 compare_effective_infall_references.py
```

This evaluator reads every HDF5 IC directly and independently reproduces the
cell-wise mode-5 rate. It does not use mode 5 to define its denominator.
For gamma=5/3 BHL tests, subsonic cases use the Bondi rate with lambda=1/4
from the asymptotic IC state; trans/supersonic cases use the
Foglizzo-Ruffert interpolation as a literature baseline rather than exact
hydrodynamic truth. Rotating and turbulent cases use the analytic continuum
IC supply (1e-4 Msun/yr) as the primary physical reference. The net resolved
10--12.5 pc shell flux of each resolution is retained separately as an IC
discretization diagnostic. After new AREPO diagnostic runs, the script also
compares logged `FFR_EFFECTIVE` against the independent IC calculation.


### Mode 5: force-corrected free-fall shell

Mode 5 remains a shell FFR estimator:

```text
Mdot_5 = C_shell * sum_i m_i / t_ff,eff,i
```

with the original central-mass free-fall time `t_ff=sqrt(r^3/(G M_central))`.
For each shell cell, AREPO's local reconstructed pressure-gradient vector is
projected along the BH-cell radial direction. The effective timescale is

```text
t_ff,eff = t_ff * sqrt(a_grav / a_net,in)
```

when `a_net,in>0`, otherwise that cell contributes zero, with

```text
a_net,in = G M_central/r^2
         + (grad P . rhat)/rho
         - v_perp^2/r
         - max(v_r,0)^2/r .
```

This uses no shell-averaged pressure profile and therefore does not assume a
spherically symmetric gas distribution. Gravity-only gas reduces exactly to
the previous shell FFR. Inward radial motion is not separately boosted, which
avoids double-counting the fact that a fast steady inflow already has less mass
resident in the shell.

Because the pressure gradient exists only after AREPO reconstructs the hydro
state, old logs cannot validate this mode. Rebuild and force-rerun:

```bash
./build_flow_tests.sh
FORCE_RERUN=1 ./run_flow_stage.sh diagnostic
python3 compare_force_ffr_references.py
```

The comparison script computes BHL references independently from the asymptotic
IC state and uses the analytic 1e-4 Msun/yr continuum supply for rotating and
turbulent ICs.
