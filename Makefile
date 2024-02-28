CC=gcc
CFLAGS=-std=c99
CFLAGS_DBG=-O0 -g $(CFLAGS)
CFLAGS_REL=-O3 -DNDEBUG $(CFLAGS)
LDFLAGS=-lm
SRCDIR=src
SRCS=$(addprefix ${SRCDIR}/,main.c)
BIN=bin
PROG=cluster

all: debug

debug: ${BIN}/${PROG}.debug
release: ${BIN}/${PROG}.release

${BIN}/${PROG}.debug: ${SRCS}
	$(CC) $(CFLAGS_DGB) $^ $(LDFLAGS) -o $@
${BIN}/${PROG}.release: ${SRCS}
	$(CC) $(CFLAGS_REL) $^ $(LDFLAGS) -o $@

.PHONY: all debug release clean
clean:
	$(RM) cluster
