all: simarrange

CFLAGS=-O3 -I . -fno-builtin-strlen  `pkg-config opencv --cflags` -DPARALLEL 
LDFLAGS=-lm `pkg-config opencv --libs` -largtable2 -ladmesh -fopenmp 

CC=gcc

simarrange: simarrange.o
simarrange.o: simarrange.c utlist.h

clean:
	rm -f simarrange simarrange.o

prefix=/usr/local

install: simarrange simarrange.1
	install -m 0755 simarrange $(prefix)/bin
	install -m 0644 simarrange.1 $(prefix)/share/man/man1

.PHONY: install
