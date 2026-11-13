bin_PROGRAMS = omf5 # omf5_dbg omf5_dbg_full omf5_nbr_data omf5_nbr_data_dbg omf5_nbr_data_dbg_full
# lib_LTLIBRARIES = libomf5.la

# DEPS_CPPFLAGS and DEPS_LDFLAGS are set by configure, flint first,
# and with antic only for flint 2.x, so that the build compiles against
# the flint headers that configure checked. The library is not checked:
# directories given in the user's LDFLAGS are searched before these.
ALL_CFLAGS = -I$(top_srcdir)/include @DEPS_CPPFLAGS@ -Wall -Werror

# The code defines its own bool (typedefs.h), which C23 makes a keyword.
# Recent autoconf selects -std=gnu23 by default, so pin the standard here.
ALL_CFLAGS += -std=gnu17

OPT_CFLAGS = -O3 -DNDEBUG
DBG_CFLAGS = -g -DDEBUG -fsanitize=address
OPT_CFLAGS += -fomit-frame-pointer

# Lines inside automake conditionals must not start with a tab,
# otherwise automake does not treat them as variable assignments.
if X86
OPT_CFLAGS += -msse3 -mavx
endif

if LINUX
DBG_CFLAGS += -pg
endif

OPT_CFLAGS += $(ALL_CFLAGS)
DBG_CFLAGS += ${ALL_CFLAGS}

DBG_FULL_CFLAGS = -DDEBUG_LEVEL_FULL
DBG_FULL_CFLAGS += ${DBG_CFLAGS}

ALL_LD_FLAGS = @DEPS_LDFLAGS@ -L${CARAT_PATH}/lib -L${CARAT_PATH}

# The libraries must come after the objects on the link line (LDADD, not LDFLAGS),
# otherwise GNU ld drops the static libfunctions.a of CARAT.
ALL_LIBS = -lfunctions @ANTIC_LIBS@ -lflint -lgmp -lm

DBG_LD_FLAGS = ${ALL_LD_FLAGS}
# DBG_LD_FLAGS += -fsanitize=leak

SOURCES = arith.c
SOURCES += aut_grp.c
SOURCES += decomposition.c
SOURCES += eigenvalues.c
SOURCES += fmpq_poly.c
SOURCES += fq_nmod_mat.c
SOURCES += fq_nmod_mpoly.c
SOURCES += fq_nmod_mpoly_mat.c
SOURCES += fq_nmod_quad.c
SOURCES += genus.c
SOURCES += hash.c
SOURCES += hecke.c
SOURCES += io.c
SOURCES += isometry.c
SOURCES += jordan.c
SOURCES += mass.c
SOURCES += matrix_tools.c
SOURCES += nbr_data.c
SOURCES += neighbor.c
SOURCES += nf_mat.c
SOURCES += orthogonalize.c
SOURCES += pivot_data.c
SOURCES += qf_inv.c
SOURCES += spinor.c
SOURCES += square_matrix.c
SOURCES += tests.c
SOURCES += weight.c


# omf5_dbg_full_SOURCES = omf5.c $(SOURCES)
# omf5_dbg_full_CFLAGS = $(DBG_FULL_CFLAGS)

# omf5_dbg_SOURCES = omf5.c $(SOURCES)
# omf5_dbg_CFLAGS = $(DBG_CFLAGS)

omf5_SOURCES = omf5.c $(SOURCES)
omf5_CFLAGS = $(OPT_CFLAGS)

# omf5_nbr_data_dbg_full_SOURCES = omf5.c $(SOURCES)
# omf5_nbr_data_dbg_full_CFLAGS = $(DBG_FULL_CFLAGS) -DNBR_DATA

# omf5_nbr_data_dbg_SOURCES = omf5.c $(SOURCES)
# omf5_nbr_data_dbg_CFLAGS = $(DBG_CFLAGS) -DNBR_DATA

# omf5_nbr_data_SOURCES = omf5.c $(SOURCES)
# omf5_nbr_data_CFLAGS = $(OPT_CFLAGS) -DNBR_DATA

# libomf5_la_SOURCES = $(SOURCES)

# libomf5_la_CFLAGS = -static $(OPT_CFLAGS)

# omf5_dbg_full_LDFLAGS = ${DBG_LD_FLAGS}

# omf5_dbg_LDFLAGS = ${DBG_LD_FLAGS}

# omf5_nbr_data_dbg_full_LDFLAGS = ${DBG_LD_FLAGS}

# omf5_nbr_data_dbg_LDFLAGS = ${DBG_LD_FLAGS}

OPT_LD_FLAGS = ${ALL_LD_FLAGS}
OPT_LD_FLAGS += -O3

omf5_LDFLAGS = ${OPT_LD_FLAGS}
omf5_LDADD = ${ALL_LIBS}
# libomf5_la_LDFLAGS = ${OPT_LD_FLAGS}

# omf5_nbr_data_LDFLAGS = ${OPT_LD_FLAGS}
