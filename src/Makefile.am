bin_PROGRAMS = A61 A61_dbg A61_dbg_full
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
SOURCES += genus.c
SOURCES += hash.c
SOURCES += hecke.c
SOURCES += jordan.c
SOURCES += mass.c
SOURCES += matrix_tools.c
SOURCES += neighbor.c
SOURCES += nf_mat.c
SOURCES += orthogonalize.c
SOURCES += tests.c
SOURCES += qf_inv.c

A61_dbg_full_SOURCES = A61.c $(SOURCES)
A61_dbg_full_CFLAGS = $(DBG_FULL_CFLAGS)

A61_dbg_SOURCES = A61.c $(SOURCES)
A61_dbg_CFLAGS = $(DBG_CFLAGS)

A61_SOURCES = A61.c $(SOURCES)
A61_CFLAGS = $(OPT_CFLAGS)

libomf5_la_SOURCES = $(SOURCES)
libomf5_la_CFLAGS = -shared -fPIC $(OPT_CFLAGS)

A61_dbg_full_LDFLAGS = ${DBG_LD_FLAGS}

A61_dbg_LDFLAGS = ${DBG_LD_FLAGS}

OPT_LDFLAGS = ${ALL_LDFLAGS}
OPT_LDFLAGS += -O3 -fPIC

A61_LDFLAGS = ${OPT_LDFLAGS}
libomf5_la_LDFLAGS = ${OPT_LDFLAGS}


