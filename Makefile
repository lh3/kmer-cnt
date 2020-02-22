CFLAGS=-g -Wall -O2
LIBS=-lz

kc-c1:kc-c1.c khashl.h ketopt.h kseq.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
