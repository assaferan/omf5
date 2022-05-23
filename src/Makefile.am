bin_PROGRAMS = A61 A61_dbg
lib_LTLIBRARIES = libomf5.la

INCLUDE_FLAGS = -I ../include

OPT_CFLAGS = -O3 -DNDEBUG -Wall -Werror -msse3 -mavx
DBG_CFLAGS = -g -Wall -Werror -DDEBUG -DDEBUG_LEVEL_FULL -fsanitize=address
OPT_CFLAGS += -fomit-frame-pointer

OPT_CFLAGS += ${INCLUDE_FLAGS}
DBG_CFLAGS += ${INCLUDE_FLAGS}

# These are not supported on macOS

if LINUX
	AM_MSG_ERROR("Still beleives this is linux...")
	OPT_CFLAGS += -fvar-tracking-assignments-toggle
	OPT_CFLAGS += -funroll-all-loops
        DBG_CFLAGS += -pg
endif

DBG_LDFLAGS = -L../libs -lfunctions -lgmp

SOURCES = arith.c
SOURCES += genus.c
SOURCES += hash.c
SOURCES += hecke.c
SOURCES += matrix_tools.c
SOURCES += neighbor.c

A61_dbg_SOURCES = A61.c $(SOURCES)
A61_dbg_CFLAGS = $(DBG_CFLAGS)

A61_SOURCES = A61.c $(SOURCES)
A61_CFLAGS = $(OPT_CFLAGS)

libomf5_la_SOURCES = $(SOURCES)
libomf5_la_CFLAGS = -shared -fPIC $(OPT_CFLAGS)

A61_dbg_LDFLAGS = ${DBG_LDFLAGS}

OPT_LDFLAGS = ${DBG_LDFLAGS}
OPT_LDFLAGS += -O3 -fPIC

A61_LDFLAGS = ${OPT_LDFLAGS}
libomf5_la_LDFLAGS = ${OPT_LDFLAGS}


