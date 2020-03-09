CFLAGS=-g -Wall -O2
CXXFLAGS=$(CFLAGS) -std=c++11
LIBS=-lz
PROG=kc-c1 kc-c2 kc-c3 kc-c4 kc-cpp1 kc-cpp2 yak-count

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

kc-c4:kc-c4.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ kc-c4.c kthread.c $(LIBS) -lpthread

yak-count:yak-count.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ yak-count.c kthread.c $(LIBS) -lpthread

kc-cpp1:kc-cpp1.cpp ketopt.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

kc-cpp2:kc-cpp2.cpp ketopt.h robin_hood.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

clean:
	rm -fr *.dSYM $(PROG)
