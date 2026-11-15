#!/bin/sh
#
# Runs the tests of omf5 (omf5 -tests): the unit tests, and the eigenforms
# and the traces of their Hecke eigenvalues for two genera.
#
# usage: scripts/run_tests.sh [path to omf5 executable]

OMF5=${1:-src/omf5}

# The eigenvalues printed to stdout are not needed here.
"$OMF5" -tests > /dev/null
rc=$?

# The exit status of omf5 -tests is not 0 when a test failed, and depends on
# the tests which failed. make check reads the status 77 as a test which was
# skipped, so every failure is returned as 1.
if [ "$rc" -ne 0 ]; then
    printf '%s -tests failed with the exit status %s\n' "$OMF5" "$rc" >&2
    exit 1
fi
