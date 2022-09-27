VPATH=src
CC=clang
CFLAGS=-I${GMP_PATH}/include -I${ANTIC_PATH}/include -Iinclude
LDFLAGS=-L${GMP_PATH}/lib -L${ANTIC_PATH}/lib -lgmp -lflint -lantic -lfunctions

TARGET_DIR = bin
TARGET = ${TARGET_DIR}/omf5

SOURCES = $(wildcard $(VPATH)/*.c)
MODULES = $(patsubst $(VPATH)/%.c,%,$(SOURCES))
OBJECTS = $(patsubst $(VPATH)/%.c,$(TARGET_DIR)/%.o,$(SOURCES))

${TARGET_DIR}/%.o : %.c
	$(CC) $(CFLAGS) -c $(patsubst $(TARGET_DIR)/%.o, $(VPATH)/%.c, $@) -o $@

${TARGET} : $(OBJECTS)
	 $(CC) $(LDFLAGS) $(OBJECTS) -o ${TARGET}

clean :
	rm -f ${TARGET_DIR}/*.o ${TARGET}
