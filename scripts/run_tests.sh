#!/bin/sh
#
# Runs the standard tests of omf5 (omf5 -tests), which check the
# eigenforms and the traces of their Hecke eigenvalues for two genera.
#
# usage: scripts/run_tests.sh [path to omf5 executable]

OMF5=${1:-src/omf5}

# The eigenvalues printed to stdout are not needed here.
"$OMF5" -tests > /dev/null
