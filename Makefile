CFLAGS=-g -Wall -O2
LIBS=-lz
PROG=kc-c1 kc-c2 kc-c3

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.PHONY:all clean

all:$(PROG)

kc-c1:kc-c1.c khashl.h ketopt.h kseq.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

kc-c2:kc-c2.c khashl.h ketopt.h kseq.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

kc-c3:kc-c3.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ kc-c3.c kthread.c $(LIBS) -lpthread

clean:
	rm -fr *.dSYM $(PROG)
