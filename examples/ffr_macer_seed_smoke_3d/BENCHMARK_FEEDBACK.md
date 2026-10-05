# Convergence-benchmark feedback models

The feedback selector is independent of the resolved accretion estimator:

- `BHBenchmarkFeedbackModel = 0`: no feedback.
- `BHBenchmarkFeedbackModel = 1`: IllustrisTNG thermal/kinetic feedback law.
- `BHBenchmarkFeedbackModel = 2`: validated FFR-MACER mechanical wind+jet feedback.

The supported backend combinations are intentionally strict:

- TNG feedback requires `BHBenchmarkAccretionTarget = 1` (direct).
- MACER feedback requires `BHBenchmarkAccretionTarget = 0` (reservoir).
- `NONE` may be used with either target for isolated tests.

This produces the six headline configurations without duplicating the
accretion implementation:

| Configuration | Accretion model | Target | Feedback |
| --- | --- | --- | --- |
| TNG_BONDI | 0 | direct | TNG |
| TNG_BOOSTED | 1 | direct | TNG |
| TNG_AMBONDI | 2 | direct | TNG |
| TNG_FFR | 3 | direct | TNG |
| MACER_BONDI | 0 | reservoir | MACER |
| FFR_MACER | 3 | reservoir | MACER |

## TNG mode switch

For the TNG feedback law,

[
chi = minleft[chi_0
  left(rac{M_{m BH}}{10^8,M_odot}ight)^eta,
  chi_{max}ight].
]

The fiducial parameters are

[
chi_0=0.002,qquad eta=2,qquad chi_{max}=0.1.
]

The feedback state is selected with the **uncapped** benchmark estimator,

[
f_{m Edd,raw} =
rac{dot M_{m raw}}{dot M_{m Edd}}.
]

The thermal state is used when (f_{m Edd,raw}gechi); otherwise the
kinetic state is used.  For model 0 this reproduces the TNG choice based on
uncapped Bondi accretion.  Models 1--3 use the same rule applied to their own
raw estimator, which is the controlled cross-model generalization.

## TNG energy generation

The energy source uses the benchmark operational accretion rate (after the
direct-backend Eddington cap), not the finite-step realized gas sink:

[
dot E_{m thermal} =
epsilon_{m f,high}epsilon_rdot M_{m operational}c^2,
]

with fiducial

[
epsilon_{m f,high}=0.1,qquad epsilon_r=0.2.
]

The kinetic state uses

[
dot E_{m kinetic} =
epsilon_{m f,kin}dot M_{m operational}c^2,
]

where

[
epsilon_{m f,kin}
=
minleft[
rac{n_{m H}}
     {f_{m thresh} n_{m H,SF}},
epsilon_{m f,kin,max}
ight].
]

Fiducial parameters are

[
f_{m thresh}=0.05,qquad
n_{m H,SF}=0.1,{m cm^{-3}},qquad
epsilon_{m f,kin,max}=0.2.
]

The density is measured with a cubic-spline kernel over the controlled fixed
proper `BHFeedbackRadius`.

## TNG kinetic bursts

Kinetic energy accumulates persistently until

[
E_{m kin,buf}ge E_{m inj,min},
qquad
E_{m inj,min}=
f_{m re}rac12 M_{m enc}sigma_{m DM}^2,
]

with fiducial (f_{m re}=20), full gas mass (M_{m enc}) inside the
feedback aperture, and the cached one-dimensional local dark-matter velocity
dispersion.

When an event fires, all available kinetic-buffer energy is released.  A
single reproducible pseudorandom unit vector is chosen for the event.  Each
active target cell receives

[
Delta {f p}_j =
m_jsqrt{rac{2Delta E,W_j}{ho_{m active}}},hat{f n}.
]

The direction changes between events.  This deliberately does **not** impose
zero momentum per event: that is the source TNG prescription.  The quadratic
kinetic term is checked to sum to (Delta E); the actual lab-frame total
energy change also contains the pre-existing velocity cross term and is
logged separately.

## Controlled aperture and active-cell policy

The original TNG cosmological implementation adapts its BH smoothing region to
a target neighbour count.  The headline convergence benchmark instead keeps
the same fixed proper feedback aperture for TNG and MACER.  This is
intentional: it isolates the feedback prescription from neighbourhood-size
changes.  Therefore the paper should call this the **TNG feedback law in the
controlled fixed aperture**, not a byte-for-byte native TNG implementation.

Public AREPO also lacks the private TNG wake-up machinery.  Only synchronized
active hydro cells are modified.  The full aperture is still used read-only
for (M_{m enc}), density, and the kinetic threshold.  If the active subset
is too small, the energy remains in its persistent buffer.  When injection is
allowed, kernel weights are renormalized over active targets so buffered
energy is not lost.

## MACER path

`BHBenchmarkFeedbackModel=2` calls the existing validated FFR-MACER feedback
implementation without changing its physical closure:

- broad bipolar wind around the coherent disc axis,
- narrow bipolar jet around the persistent jet axis,
- separate wind and jet energy reservoirs,
- wind rest-mass return,
- binding/virial burst thresholds,
- active-target/lobe safeguards,
- exact zero-net-momentum bipolar packets,
- exact requested kinetic-energy increment,
- packet backlog and feedback-aware timestep handling.

Thus `FFR + reservoir + MACER` remains the validated production model, while
`Bondi + reservoir + MACER` changes only the resolved supply estimator.

## Snapshot diagnostics

The benchmark writes these additional Type-5 fields:

- `BH_BenchmarkMdotRaw`
- `BH_BenchmarkMdotOperational`
- `BH_TNGMode`
- `BH_TNGRawEddRatio`
- `BH_TNGModeThreshold`
- `BH_TNGFeedbackPower`
- `BH_TNGThermalEnergyBuffer`
- `BH_TNGKineticEnergyBuffer`
- `BH_TNGKineticThresholdEnergy`

The existing MACER diagnostics remain unchanged.

## Restart format

Adding persistent TNG energy buffers changes the benchmark BHP layout.  Native
benchmark restarts therefore use FFR restart version 3 and explicitly reject
older incompatible native restart records.  This does not modify the
validated `ffr-macer` branch or its restart format.
