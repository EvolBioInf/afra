CFLAGS= -W -Wall -g -O0 -fsanitize=address -fms-extensions
CPPFLAGS= -std=gnu99

.PHONY: all clean check

all: alfsv

alfsv: alfsv.o io.o matrix.o graph.o
	$(CC) $(CFLAGS) -o $@ $^

check: alfsv
	./alfsv ../../data/eco29.mat

clean:
	rm -f alfsv *.o
