CC=gcc
CFLAGS=-std=gnu99 -Wall -Wextra
CFLAGS_DBG=-O0 -g $(CFLAGS)
CFLAGS_REL=-O3 -DNDEBUG -g $(CFLAGS)
CFLAGS_PROF=-O3 -DNDEBUG -pg -g $(CFLAGS)
LDFLAGS=-lm
SRCDIR=src
SRCS=$(addprefix ${SRCDIR}/,main.c)
BIN=bin
PROG=cluster

NROWS?=2
NCOLS?=2
NT?=64
CFLAGS_CFG=-DNROWS=$(NROWS) -DNCOLS=$(NCOLS) -DNT=$(NT)

all: release

debug: ${BIN}/${PROG}.$(NT)_$(NROWS)_$(NCOLS).debug
release: ${BIN}/${PROG}.$(NT)_$(NROWS)_$(NCOLS).release
prof: ${BIN}/${PROG}.$(NT)_$(NROWS)_$(NCOLS).prof

${BIN}/${PROG}.$(NT)_$(NROWS)_$(NCOLS).debug: ${SRCS}
	$(CC) $(CFLAGS_DGB) $(CFLAGS_CFG) $^ $(LDFLAGS) -o $@
${BIN}/${PROG}.$(NT)_$(NROWS)_$(NCOLS).release: ${SRCS}
	$(CC) $(CFLAGS_REL) $(CFLAGS_CFG) $^ $(LDFLAGS) -o $@
${BIN}/${PROG}.$(NT)_$(NROWS)_$(NCOLS).prof: ${SRCS}
	$(CC) $(CFLAGS_PROF) $(CFLAGS_CFG) $^ $(LDFLAGS) -o $@


.PHONY: all debug release prof clean
clean:
	$(RM) ${BIN}/*
