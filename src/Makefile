CFLAGS= -W -Wall -g -Os -fsanitize=address -fms-extensions
CPPFLAGS= -std=gnu99

.PHONY: all clean

all: alfsv


alfsv: alfsv.o io.o matrix.o graph.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f alfsv *.o
