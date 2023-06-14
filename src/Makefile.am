bin_PROGRAMS = omf5 # omf5_dbg omf5_dbg_full omf5_nbr_data omf5_nbr_data_dbg omf5_nbr_data_dbg_full
# lib_LTLIBRARIES = libomf5.la

ALL_CFLAGS = -I ../include -I${GMP_PATH}/include -I${ANTIC_PATH}/include -I${FLINT_PATH}/include -Wall -Werror

OPT_CFLAGS = -O3 -DNDEBUG -msse3 -mavx
DBG_CFLAGS = -g -DDEBUG -fsanitize=address
OPT_CFLAGS += -fomit-frame-pointer

if LINUX
	OPT_CFLAGS += -fvar-tracking-assignments-toggle
	OPT_CFLAGS += -funroll-all-loops
        DBG_CFLAGS += -pg
endif

if OSX
	ALL_CFLAGS += -mmacosx-version-min=10.10
endif

OPT_CFLAGS += $(ALL_CFLAGS)
DBG_CFLAGS += ${ALL_CFLAGS}

DBG_FULL_CFLAGS = -DDEBUG_LEVEL_FULL
DBG_FULL_CFLAGS += ${DBG_CFLAGS}

ALL_LDFLAGS = -L${GMP_PATH}/lib -L${ANTIC_PATH}/lib -L${FLINT_PATH}/lib -lfunctions -lgmp -lantic -lflint -lm

DBG_LD_FLAGS = ${ALL_LDFLAGS}
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

OPT_LDFLAGS = ${ALL_LDFLAGS}
OPT_LDFLAGS += -O3

omf5_LDFLAGS = ${OPT_LDFLAGS}
# libomf5_la_LDFLAGS = ${OPT_LDFLAGS}

# omf5_nbr_data_LDFLAGS = ${OPT_LDFLAGS}
