VPATH = src
CC = clang
CFLAGS = -I${GMP_PATH}/include -I${ANTIC_PATH}/include -Iinclude
LDFLAGS = -L${GMP_PATH}/lib -L${ANTIC_PATH}/lib -lgmp -lflint -lantic -lfunctions

DBG_CFLAGS = -g -DDEBUG #-fsanitize=address

DBG_FULL_CFLAGS = -DDEBUG_LEVEL_FULL $(DBG_CFLAGS)

TARGET_DIR = bin
TARGET = ${TARGET_DIR}/omf5

SOURCES = $(wildcard $(VPATH)/*.c)
MODULES = $(patsubst $(VPATH)/%.c,%,$(SOURCES))
OBJECTS = $(patsubst $(VPATH)/%.c,$(TARGET_DIR)/%.o,$(SOURCES))

DBG_OBJECTS = $(patsubst %.o, %_dbg.o, $(OBJECTS))
DBG_TARGET = ${TARGET}_dbg

DBG_FULL_OBJECTS = $(patsubst %.o, %_dbg_full.o, $(OBJECTS))
DBG_FULL_TARGET = ${TARGET}_dbg_full

all: ${TARGET} ${DBG_TARGET} ${DBG_FULL_TARGET}

${TARGET_DIR}/%.o : %.c
	$(CC) $(CFLAGS) -c $(patsubst $(TARGET_DIR)/%.o, $(VPATH)/%.c, $@) -o $@

${TARGET_DIR}/%_dbg.o : %.c
	$(CC) $(CFLAGS) $(DBG_CFLAGS) -c $(patsubst $(TARGET_DIR)/%_dbg.o, $(VPATH)/%.c, $@) -o $@

${TARGET_DIR}/%_dbg_full.o : %.c
	$(CC) $(CFLAGS) $(DBG_FULL_CFLAGS) -c $(patsubst $(TARGET_DIR)/%_dbg_full.o, $(VPATH)/%.c, $@) -o $@

${TARGET} : $(OBJECTS)
	 $(CC) $(LDFLAGS) $(OBJECTS) -o ${TARGET}

${DBG_TARGET}: $(DBG_OBJECTS)
	 $(CC) $(LDFLAGS) $(DBG_OBJECTS) -o ${DBG_TARGET}

${DBG_FULL_TARGET}: $(DBG_FULL_OBJECTS)
	 $(CC) $(LDFLAGS) $(DBG_FULL_OBJECTS) -o ${DBG_FULL_TARGET}

clean :
	rm -f ${TARGET_DIR}/*.o ${TARGET}
	rm -f ${TARGET_DIR}/*.o_dbg ${DBG_TARGET}
	rm -f ${TARGET_DIR}/*.o_dbg ${DBG_FULL_TARGET}
