SUBDIRS = src

ACLOCAL_AMFLAGS = -I m4

# make check runs omf5 -tests and checks the README examples.
TESTS = scripts/run_tests.sh scripts/check_examples.sh
EXTRA_DIST = $(TESTS)

# The headers are not listed in src/Makefile.am, so distribute them here.
EXTRA_DIST += include
