SUBDIRS = src

ACLOCAL_AMFLAGS = -I m4

# make check runs omf5 -tests and checks the README examples.
TESTS = scripts/run_tests.sh scripts/check_examples.sh
EXTRA_DIST = $(TESTS)

# The headers are not listed in src/Makefile.am, so distribute them here.
EXTRA_DIST += include

# automake only distributes the license on its own if it is named COPYING.
EXTRA_DIST += LICENSE.txt

# EXTRA_DIST copies the whole include directory, so remove anything that is
# not a header, such as .DS_Store or editor backups.
dist-hook:
	find $(distdir)/include -type f ! -name '*.h' -exec rm -f {} +
