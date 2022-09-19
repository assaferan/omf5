bin_PROGRAMS = omf5 omf5_dbg omf5_dbg_full
lib_LTLIBRARIES = libomf5.la

ALL_CFLAGS = -I ../include -Wall -Werror

OPT_CFLAGS = -O3 -DNDEBUG -msse3 -mavx
DBG_CFLAGS = -g -DDEBUG -fsanitize=address
OPT_CFLAGS += -fomit-frame-pointer

# These are not supported on macOS

if LINUX
	AM_MSG_ERROR("Still beleives this is linux...")
	OPT_CFLAGS += -fvar-tracking-assignments-toggle
	OPT_CFLAGS += -funroll-all-loops
        DBG_CFLAGS += -pg
endif

if OSX
	ALL_CFLAGS += -mmacosx-version-min=10.10
endif

OPT_CFLAGS += ${ALL_CFLAGS}
DBG_CFLAGS += ${ALL_CFLAGS}

DBG_FULL_CFLAGS = -DDEBUG_LEVEL_FULL
DBG_FULL_CFLAGS += ${DBG_CFLAGS}

ALL_LDFLAGS = -lfunctions -lgmp -lantic -lflint

DBG_LD_FLAGS = ${ALL_LDFLAGS}
# DBG_LD_FLAGS += -fsanitize=leak

SOURCES = arith.c
SOURCES += fq_nmod_mat.c
SOURCES += fq_nmod_quad.c
SOURCES += fq_nmod_poly_mat.c
SOURCES += genus.c
SOURCES += hash.c
SOURCES += hecke.c
SOURCES += jordan.c
SOURCES += mass.c
SOURCES += matrix_tools.c
SOURCES += nbr_data.c
SOURCES += neighbor.c
SOURCES += nf_mat.c
SOURCES += orthogonalize.c
SOURCES += pivot_data.c
SOURCES += tests.c
SOURCES += qf_inv.c

omf5_dbg_full_SOURCES = omf5.c $(SOURCES)
omf5_dbg_full_CFLAGS = $(DBG_FULL_CFLAGS)

omf5_dbg_SOURCES = omf5.c $(SOURCES)
omf5_dbg_CFLAGS = $(DBG_CFLAGS)

omf5_SOURCES = omf5.c $(SOURCES)
omf5_CFLAGS = $(OPT_CFLAGS)

libomf5_la_SOURCES = $(SOURCES)
libomf5_la_CFLAGS = -shared -fPIC $(OPT_CFLAGS)

omf5_dbg_full_LDFLAGS = ${DBG_LD_FLAGS}

omf5_dbg_LDFLAGS = ${DBG_LD_FLAGS}

OPT_LDFLAGS = ${ALL_LDFLAGS}
OPT_LDFLAGS += -O3 -fPIC

omf5_LDFLAGS = ${OPT_LDFLAGS}
libomf5_la_LDFLAGS = ${OPT_LDFLAGS}


