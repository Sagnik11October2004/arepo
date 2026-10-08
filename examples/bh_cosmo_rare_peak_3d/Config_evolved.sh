#!/bin/bash

# Post-split cosmological executable for z<=22 snapshots that already contain
# gas.  It deliberately omits GENERATE_GAS_IN_ICS.

SELFGRAVITY
HIERARCHICAL_GRAVITY
CELL_CENTER_GRAVITY
ALLOW_DIRECT_SUMMATION
DIRECT_SUMMATION_THRESHOLD=500
PMGRID=64
RCUT=5.0

NSOFTTYPES=2
MULTIPLE_NODE_SOFTENING
ADAPTIVE_HYDRO_SOFTENING
TREE_BASED_TIMESTEPS

REGULARIZE_MESH_CM_DRIFT
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED

DOUBLEPRECISION=1
DOUBLEPRECISION_FFTW
INPUT_IN_DOUBLEPRECISION
OUTPUT_IN_DOUBLEPRECISION
OUTPUT_COORDINATES_IN_DOUBLEPRECISION
NGB_TREE_DOUBLEPRECISION

HAVE_HDF5
COOLING

BLACKHOLE_FFR
# Controlled post-seed native-restart experiment (1 physical Myr dormancy).
# Reference epoch: the common post-seed checkpoint at z ~= 20.88374.
BH_FFR_DELAY_QUIET_MYR=1.0
BH_FFR_DELAY_QUIET_START_Z=20.88374
BH_FFR_SEED_ONLY_MOST_MASSIVE
BH_FFR_SEED_MAX_REDSHIFT=22

FOF
FOF_PRIMARY_LINK_TYPES=2
FOF_SECONDARY_LINK_TYPES=1+32

PROCESS_TIMES_OF_OUTPUTLIST
REDUCE_FLUSH

# Science-branch policy:
# - shell FFR production comparisons use A_ff=1e-2 and 1e-1 at runtime;
# - mode 5 (ConvJ/Mach shell FFR) uses its unit geometric-shell base with
#   resolved f_M*f_j suppression and no extra 1e-3 prefactor;
# - timestep floors remain the standard AREPO runtime MinSizeTimestep only.
#   No BLACKHOLE_FFR-specific minimum-timestep floor is compiled here.
