#!/bin/make -f

CC=clang

CFLAGS= -W -Wall -g -O1 -fsanitize=address -fms-extensions -Wno-microsoft -Wno-missing-field-initializers
CPPFLAGS= -std=c11

.PHONY: all clean check

all: afra

afra: afra.o io.o matrix.o graph.o quartet.o consense.o
	$(CC) $(CFLAGS) -o $@ $^

check: afra
	./afra ../../data/eco29.mat

clean:
	rm -f afra *.o
