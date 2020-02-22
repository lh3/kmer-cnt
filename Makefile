CFLAGS=-g -Wall -O2
LIBS=-lz
PROG=kc-c1

.PHONY:clean

kc-c1:kc-c1.c khashl.h ketopt.h kseq.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

clean:
	rm -fr *.dSYM $(PROG)
