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
