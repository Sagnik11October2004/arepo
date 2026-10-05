# Convergence-benchmark feedback models

The feedback selector is independent of the resolved accretion estimator:

- `BHBenchmarkFeedbackModel = 0`: no feedback.
- `BHBenchmarkFeedbackModel = 1`: IllustrisTNG thermal/kinetic feedback law.
- `BHBenchmarkFeedbackModel = 2`: validated FFR-MACER mechanical wind+jet feedback.

The allowed backend combinations are intentionally strict:

- TNG feedback requires `BHBenchmarkAccretionTarget = 1` (direct).
- MACER feedback requires `BHBenchmarkAccretionTarget = 0` (reservoir).
- `NONE` may be used with either target for isolated tests.

This yields the six headline configurations without duplicating the accretion
implementation:

| Configuration | Accretion model | Target | Feedback |
| --- | --- | --- | --- |
| TNG_BONDI | 0 | direct | TNG |
| TNG_BOOSTED | 1 | direct | TNG |
| TNG_AMBONDI | 2 | direct | TNG |
| TNG_FFR | 3 | direct | TNG |
| MACER_BONDI | 0 | reservoir | MACER |
| FFR_MACER | 3 | reservoir | MACER |

## TNG mode switch

The TNG feedback switch is

\[
\chi =
\min\left[
\chi_0
\left(\frac{M_{\rm BH}}{10^8\,M_\odot}\right)^\beta,
\chi_{\max}
\right],
\]

with fiducial values

\[
\chi_0=0.002,\qquad \beta=2,\qquad \chi_{\max}=0.1.
\]

The feedback state is selected from the uncapped benchmark estimator,

\[
f_{\rm Edd,raw} =
\frac{\dot M_{\rm raw}}{\dot M_{\rm Edd}}.
\]

Thermal mode is used for \(f_{\rm Edd,raw}\ge\chi\), otherwise kinetic mode is
used. For model 0 this is the native TNG choice based on uncapped Bondi
accretion. Models 1--3 use the same switch applied to their own raw estimator;
that extension is part of the controlled cross-model benchmark.

## TNG energy generation

The energy source uses the benchmark operational direct accretion rate, i.e.
the selected estimator after the direct-backend Eddington cap. It is kept
separate from the finite-step realized conservative gas sink.

Thermal mode:

\[
\dot E_{\rm thermal}
=
\epsilon_{\rm f,high}\epsilon_r
\dot M_{\rm operational}c^2,
\]

with

\[
\epsilon_{\rm f,high}=0.1,\qquad \epsilon_r=0.2.
\]

Kinetic mode:

\[
\dot E_{\rm kinetic}
=
\epsilon_{\rm f,kin}
\dot M_{\rm operational}c^2,
\]

where

\[
\epsilon_{\rm f,kin}
=
\min\left[
\frac{n_{\rm H}}
     {f_{\rm thresh}n_{\rm H,SF}},
\epsilon_{\rm f,kin,max}
\right].
\]

The fiducial benchmark values are

\[
f_{\rm thresh}=0.05,\qquad
n_{\rm H,SF}=0.1065\,{\rm cm^{-3}},\qquad
\epsilon_{\rm f,kin,max}=0.2.
\]

The star-formation density is a configurable physical-\(n_{\rm H}\) proxy in
this public-AREPO benchmark. The checked-in value, \(0.1065\,{\rm cm^{-3}}\),
is the physical TNG threshold corresponding to the commonly quoted rounded
value \(\simeq0.1\,{\rm cm^{-3}}\). It should be matched to the star-formation
model if a later production setup deliberately changes that threshold.

## TNG kinetic bursts

Kinetic energy is stored persistently until

\[
E_{\rm kin,buf}\ge E_{\rm inj,min},
\qquad
E_{\rm inj,min}
=
f_{\rm re}\frac{1}{2}M_{\rm enc}\sigma_{\rm DM}^2,
\]

with the Weinberger et al. fiducial \(f_{\rm re}=20\). \(M_{\rm enc}\) is the
full gas mass inside the feedback aperture and \(\sigma_{\rm DM}\) is the
cached one-dimensional local dark-matter velocity dispersion.

When an event fires, the available kinetic buffer is released. A single
reproducible pseudorandom unit vector is chosen for that event and changes
between events. Each active target receives

\[
\Delta{\bf p}_j
=
m_j
\sqrt{
\frac{2\Delta E\,W_j}{\rho_{\rm active}}
}
\,\hat{\bf n}.
\]

This deliberately does not force zero net momentum for an individual event:
that is the TNG prescription. The code verifies that the sum of the quadratic
kinetic terms equals \(\Delta E\), and separately logs the actual lab-frame
energy change and injected momentum. In TNG these quantities conserve only in
the average over randomly reoriented events.

## Controlled aperture and active-cell policy

The cosmological TNG implementation normally adapts the BH neighbourhood to a
target neighbour count. The convergence benchmark instead keeps the same fixed
proper `BHFeedbackRadius` for TNG and MACER so changes in the headline
comparison come from the feedback prescription, not a different numerical
aperture. Weinberger et al. also used a fixed physical feedback radius in their
idealized kinetic-wind convergence tests.

The paper should therefore describe this implementation as the **TNG feedback
law in a controlled fixed aperture**, rather than as a byte-for-byte native TNG
BH implementation.

The public AREPO base used here does not provide the private TNG wake-up
machinery. Only synchronized active hydro cells are modified. All gas inside
the aperture is still used read-only for \(M_{\rm enc}\), the kernel density,
and the kinetic threshold. If the active subset is too small, energy remains
buffered. Injection weights are renormalized over the active subset so no
stored feedback energy is silently discarded.

## MACER path

`BHBenchmarkFeedbackModel=2` calls the already validated FFR-MACER feedback
implementation. Its physical closure is not rewritten:

- broad bipolar wind around the coherent disc axis;
- narrow bipolar jet around the persistent jet axis;
- separate wind and jet energy reservoirs;
- wind rest-mass return;
- binding/virial burst thresholds;
- active-target and lobe safeguards;
- zero-net-momentum bipolar packets;
- exact requested packet kinetic-energy increment;
- backlog and feedback-aware timestep handling.

Thus `FFR + reservoir + MACER` remains the validated production-model
combination, while `Bondi + reservoir + MACER` changes only the resolved
supply estimator.

## Snapshot diagnostics

Additional Type-5 output fields are:

- `BH_BenchmarkMdotRaw`
- `BH_BenchmarkMdotOperational`
- `BH_BenchmarkMdotRealized`
- `BH_BenchmarkSinkRateRatio`
- `BH_BenchmarkCumulativeOperationalMass`
- `BH_BenchmarkCumulativeRealizedMass`
- `BH_BenchmarkCumulativeBHMassGrowth`
- `BH_TNGMode`
- `BH_TNGRawEddRatio`
- `BH_TNGModeThreshold`
- `BH_TNGFeedbackPower`
- `BH_TNGThermalEnergyBuffer`
- `BH_TNGKineticEnergyBuffer`
- `BH_TNGKineticThresholdEnergy`
- `BH_TNGCumulativeGeneratedEnergy`
- `BH_TNGCumulativeInjectedEnergy`
- `BH_TNGActiveTargetMassFraction`
- `BH_TNGThermalBufferAgeCodeTime`
- `BH_TNGKineticBufferAgeCodeTime`
- `BH_TNGThermalBufferAgeTransactions`
- `BH_TNGKineticBufferAgeTransactions`

The cumulative TNG model-energy ledger is checked as

[
E_{m generated}
=
E_{m injected}
+
E_{m thermal,buf}
+
E_{m kinetic,buf}.
]

For kinetic events, (E_{m injected}) is the quadratic kick energy that
defines the TNG feedback reservoir; the separately logged lab-frame energy
change can differ because of the pre-existing-velocity cross term.

The active feedback mass fraction
(M_{m active,fb}/M_{m enc,fb}), current thermal/kinetic buffers, and
buffer ages in both code time and BH transactions are logged and output.
These diagnostics are intended to quantify any delayed coupling caused by
active-only public-AREPO feedback, especially when
`BHFeedbackRadius > BHAccretionRadius`.  They do not change the coupling
algorithm.

The existing MACER diagnostics are unchanged.

## Restart format

Persistent TNG feedback buffers extend the benchmark BHP state. Native restart
files written by this branch therefore use FFR restart version 4 and reject
older incompatible records explicitly. The validated `ffr-macer` branch is
untouched and retains its own existing restart format.