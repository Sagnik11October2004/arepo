#!/bin/bash            # this line only there to enable syntax highlighting in this file

#############################################################
# FFR-MACER science configuration (version-1 baseline)      #
#                                                           #
# This is intentionally separate from the tracked Config.sh #
# used by AREPOVTK. Add problem-specific physics here or    #
# in a derived config; do not turn Config.sh into the BH    #
# production configuration.                                #
#############################################################

HAVE_HDF5
SELFGRAVITY
DOUBLEPRECISION=1

# FFR-MACER Type-5 black-hole subgrid model.
BLACKHOLE_FFR
