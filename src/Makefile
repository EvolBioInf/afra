#!/bin/make -f

CC=clang

CFLAGS= -W -Wall -g -O1 -fsanitize=address -fms-extensions -Wno-microsoft -Wno-missing-field-initializers
CPPFLAGS= -std=c11

.PHONY: all clean check

all: alfsv

alfsv: alfsv.o io.o matrix.o graph.o quartet.o consense.o
	$(CC) $(CFLAGS) -o $@ $^

check: alfsv
	./alfsv ../../data/eco29.mat

clean:
	rm -f alfsv *.o
