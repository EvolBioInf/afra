CFLAGS= -W -Wall -g -O0
CPPFLAGS= -std=gnu99


all: alfsv


alfsv: alfsv.o io.o matrix.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f alfsv *.o
