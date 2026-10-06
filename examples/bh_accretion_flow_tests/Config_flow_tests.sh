#!/bin/bash

# Controlled non-cosmological BH accretion-flow benchmark.
# Type-0 gas + one pre-existing Type-5 BH are read directly from HDF5 ICs.
# FoF remains compiled because BLACKHOLE_FFR's production seeding module
# requires it at compile time; all flow-suite output lists use DumpFlag=2,
# which suppresses FoF catalogue generation for these DM-free tests.

SELFGRAVITY

PMGRID=128
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

BLACKHOLE_FFR

FOF
FOF_PRIMARY_LINK_TYPES=2
FOF_SECONDARY_LINK_TYPES=1+32

# ICs already contain gas. Never enable GENERATE_GAS_IN_ICS here.

PROCESS_TIMES_OF_OUTPUTLIST
REDUCE_FLUSH
