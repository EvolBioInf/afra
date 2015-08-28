CFLAGS= -W -Wall -g -O2
CPPFLAGS= -std=gnu99


all: alfsv


alfsv: alfsv.o io.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f alfsv *.o
