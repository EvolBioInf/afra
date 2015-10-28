#!/bin/make -f

CFLAGS= -W -Wall -g -O1 -fsanitize=address -fms-extensions -Wno-microsoft
CPPFLAGS= -std=c11

.PHONY: all clean check

all: alfsv

alfsv: alfsv.o io.o matrix.o graph.o quadruple.o
	$(CC) $(CFLAGS) -o $@ $^

check: alfsv
	./alfsv ../../data/eco29.mat

clean:
	rm -f alfsv *.o
