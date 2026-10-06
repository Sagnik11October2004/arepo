# High-resolution reference (“truth”) accretion suite

This directory builds a separate numerical reference suite for the 14 physical
flow cases in bh_accretion_flow_tests. Its purpose is to estimate the resolved
accretion rate against which the five coarse subgrid estimators can be compared.

## What counts as the reference rate?

The reference runs do **not** define the answer from Bondi, boosted Bondi,
AM-Bondi, or the production FFR coefficient. They resolve the central flow
with nested cells and place a small rapidly draining absorber at the centre.
The measured rate is the late-time slope of cumulative gas mass swallowed by
that absorber,

\[
\dot M_{\rm ref} = dM_{\rm swallowed}/dt.
\]

The absorber is implemented through the existing direct volume-capture backend
with BHFreeFallA=100, no Eddington restriction, no feedback, and a 50% maximum
removal per active cell step. With these parameters the inner sink is intended
to be much faster than the resolved supply and therefore act as an absorbing
boundary, not as the physical accretion-timescale model being tested.

The sink radius shrinks with the finest resolution:

| level | finest dx | sink radius (=4 dx) |
|---|---:|---:|
| HR4 | 0.390625 pc | 1.5625 pc |
| HR5 | 0.1953125 pc | 0.78125 pc |
| HR6 | 0.09765625 pc | 0.390625 pc |

Only the central region is refined. The refinement ladder is 40, 25, 12.5,
6.25, 3.125, and 1.5625 pc for successive split levels, keeping the total gas
cell count modest (roughly 57k, 73k, and 88k for HR4/HR5/HR6).

A case is promoted to actual_mdot_msun_yr only when the late-time HR5 -> HR6
mean sink rate changes by no more than 10%. The analyzer also records
late-time variability and direct mass fluxes through resolved shells. Failure
to converge is reported as such; HR6 is never silently declared truth.

For BHL cases the analytic BHL rate is retained only as an independent
cross-check. For rotating and turbulent cases the resolved sink rate is the
primary reference.

## Generate

From examples/bh_accretion_flow_tests:

    python3 reference_truth/generate_reference_ics.py
    python3 reference_truth/generate_reference_runs.py
    python3 reference_truth/verify_reference_suite.py

This creates 42 reference ICs and 42 runs:

- 14 physical cases,
- HR4, HR5 and HR6 for each case.

Generated binary ICs, parameter files, outputs, logs and result tables are
git-ignored.

## Run

All reference runs sequentially:

    ./reference_truth/run_reference_truth.sh

Or one selected case/level:

    FAMILY=bhl TAG=M5p0 LEVEL=HR6 ./reference_truth/run_reference_truth.sh

The default is 16 MPI ranks per reference run. EXEC and MPIEXEC may be
overridden in the environment.

Do not run this concurrently on the same laptop cores with the existing
diagnostic matrix unless oversubscription is intentional. Generating the ICs
and parameter files while other runs execute is safe.

## Analyze

After any subset has finished:

    python3 reference_truth/analyze_reference_truth.py

Outputs:

- reference_truth/reference_truth_runs.csv
- reference_truth/reference_truth_cases.csv
- reference_truth/reference_truth_summary.json

After a physical case has a converged HR5/HR6 reference and its coarse
diagnostic logs exist:

    python3 compare_models_to_truth.py

This writes models_vs_actual.csv and reports

\[
\dot M_{\rm model}/\dot M_{\rm actual}
\]

for TNG, boosted BHL, AM-Bondi, volume FFR and shell FFR at L0/L1/L2.

For FFR, the main comparison uses the configured rate including A_ff=1e-3.
Extra CSV columns also report FFR divided by A_ff * dotM_actual to diagnose
only the geometry/free-fall part of the model.

## Interpretation and remaining validation

Resolution/sink-radius convergence is the primary truth criterion. Before
using these rates in a publication-quality claim, also inspect the recorded
radial fluxes and verify that changing the absorber strength/cap in a few
representative HR6 runs does not change the mean supplied rate. If it does,
the absorber is still influencing the solution and the corresponding
actual_mdot should not be treated as physical truth.
