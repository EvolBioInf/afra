CFLAGS= -W -Wall -g -O0 -fsanitize=address
CPPFLAGS= -std=gnu99


all: alfsv


alfsv: alfsv.o io.o matrix.o graph.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f alfsv *.o
