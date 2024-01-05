CC=gcc
CFLAGS=-O0 -g -std=c99
LDFLAGS=-lm
SRCDIR=src
SRCS=$(addprefix ${SRCDIR}/,main.c)
BIN=bin

${BIN}/cluster: ${SRCS}
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

.PHONY: clean
clean:
	$(RM) cluster
