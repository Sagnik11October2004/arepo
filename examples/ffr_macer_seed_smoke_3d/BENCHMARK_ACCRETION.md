# Accretion convergence benchmark layer

This branch starts from the validated FFR-MACER implementation and adds a
runtime-selectable accretion comparison layer.  The original `ffr-macer`
branch is not modified.

## Common resolved environment

All benchmark models currently use the same fixed proper
`BHAccretionRadius`.  Gas cells whose mesh-generating points lie inside the
aperture contribute to a common top-hat estimate:

- physical density: total enclosed gas mass divided by the summed physical
  cell volume;
- sound speed: volume-weighted `get_sound_speed()`;
- bulk peculiar velocity: volume-weighted gas velocity;
- relative velocity: magnitude of bulk gas velocity minus BH velocity;
- circulation speed:
  `Vphi = |sum m_i (r_i x v_rel,i)| / (Mgas Racc)`;
- hydrogen number density:
  `nH = XH rho / mp`, with the code value `XH=0.76`.

Comoving quantities are converted before the accretion algebra:
`rho_phys = rho_comoving/a^3`, proper offsets are `a dr`, and peculiar
velocities are `P.Vel/a`.

## Accretion models

Set `BHBenchmarkAccretionModel` to:

- `0`: TNG-like unboosted Bondi,
  `Mdot = 4 pi G^2 Mbh^2 rho / cs^3`.
- `1`: boosted Bondi-Hoyle-Lyttleton,
  `Mdot = alpha 4 pi G^2 Mbh^2 rho / (cs^2+vrel^2)^(3/2)`.
- `2`: angular-momentum-regulated Bondi-Hoyle-Lyttleton,
  `Mdot = Mdot_BHL min[(cs/Vphi)^3/Cvisc,1]`.
- `3`: the validated resolved FFR estimator.

For boosted Bondi, `BHBenchmarkBoostMode=0` uses a constant
`BHBenchmarkBoostAlpha`.  Mode 1 uses the Booth-Schaye form
`alpha=1` below `BHBenchmarkBoostDensityThreshold` and
`alpha=(nH/nstar)^beta` above it.

The default angular-momentum parameter is the EAGLE reference value
`Cvisc=2*pi`.

## Backend target

Set `BHBenchmarkAccretionTarget` to:

- `0`: send captured mass to the unresolved reservoir.  No Eddington cap is
  imposed on the resolved supply.  This is the backend required for
  FFR-MACER and the future Bondi+MACER ablation.
- `1`: accrete directly into `BHMass`.  The operational rate is capped at
  `BHBenchmarkEddingtonFactor * MdotEdd`, with
  `BHBenchmarkRadiativeEfficiency` used in the Eddington rate.  The default
  values are 1 and 0.2, respectively.

Direct mode intentionally bypasses the FFR reservoir, MACER state machine,
wind and jet injection in this iteration.  TNG thermal/kinetic feedback is
not implemented yet.

## Conservative sink

All four models use the existing overlap-safe two-pass cell sink after a
common environment pass.

For FFR, the radial free-fall coefficient is retained cell by cell.  If the
direct backend Eddington-limits the raw FFR estimate, every FFR cell
coefficient is multiplied by the same rate ratio.

For Bondi-family models, the BH-level operational rate is converted to a
uniform aperture coefficient

`lambda = Mdot_operational / Mgas`.

The existing exact finite-step cell sink, per-cell sink cap, overlap
partition, gas momentum removal and global mass/momentum ledger are then
used unchanged.

The runtime log line beginning with `BH_BENCHMARK: accretion` records the
common environment, raw model rate, Eddington rate, operational rate, boost
factor and angular-momentum limiter.

## Default compatibility

The checked-in smoke-test parameter file uses

`BHBenchmarkAccretionModel = 3`
`BHBenchmarkAccretionTarget = 0`
`BHBenchmarkFeedbackModel = 2`

which recovers the validated FFR-MACER capture/reservoir path.  The added
common-environment pass changes only diagnostic/search work, not the FFR
capture coefficient or mass partition.
