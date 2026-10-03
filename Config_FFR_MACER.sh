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

# FoF is executed at full synchronization points for conservative BH seeding.
# Type-1 dark matter defines the primary FoF groups; gas (Type 0) and FFR
# black holes (Type 5) are attached to the nearest primary group.
FOF
FOF_PRIMARY_LINK_TYPES=2
FOF_SECONDARY_LINK_TYPES=1+32
