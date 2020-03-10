## Getting Started

```sh
git clone https://github.com/lh3/kmer-cnt
cd kmer-cnt
make  # C++11 required to compile the two C++ implementations
wget https://github.com/lh3/kmer-cnt/releases/download/v0.1/M_abscessus_HiSeq_10M.fa.gz
./yak-count M_abscessus_HiSeq_10M.fa.gz > kc-c4.out
```

## Introduction

K-mer counting is the foundation of many mappers, assemblers and miscellaneous
tools (e.g. genotypers, metagenomics profilers, etc). It is one of the most
important classes of algorithms in Bioinformatics. Here we will implement basic
k-mer counting algorithms but with advanced engineering tricks. We will see how
far better engineering can go.

In this repo, each `{kc,yak}-*.*` file implements a standalone k-mer counter.
As to other files: ketopt.h is a command line option parser; khashl.h is a
generic hash table library in C; kseq.h is a fasta/fastq parser; kthread.{h,c}
provides two multi-threading models; robin\_hood.h is a C++11 hash table
library.

## Results

We provide eight k-mer counters, which are detailed below the result table. All
implementations count canonical k-mers, the lexicographically *smaller* k-mer
between the k-mers on the two DNA strands.

The following table shows the timing and peak memory of different
implementations for counting 31-mers from 2.5 million pairs of 100bp reads
sampled from the HiSeq *M. abscessus* 6G-0125-R dateset in [GAGE-B][gage-b].
They were run on a Linux server equipped with two EPYC 7301 CPUs and 512GB RAM.

|Implementation                 |Limitation          |Elapsed time (s)|CPU time (s)|Peak RAM (GB)|
|:------------------------------|:-------------------|---------------:|-----------:|------------:|
|[kc-py1](kc-py1.py) + Python3.7|                    |           499.6|       499.5|         8.15|
|[kc-py1](kc-py1.py) + Pypy7.3  |                    |          1220.8|      1220.8|        12.21|
|[kc-cpp1](kc-cpp1.cpp)         |                    |           528.0|       527.9|         8.27|
|[kc-cpp2](kc-cpp2.cpp)         |                    |           319.6|       319.6|         6.90|
|[kc-c1](kc-c1.c)               |<=32-mer            |            39.3|        38.3|         1.52|
|[kc-c2](kc-c2.c)               |<=32-mer; <1024 count|           38.7|        37.9|         1.05|
|[kc-c3](kc-c3.c)               |<=32-mer; <1024 count|           34.1|        38.7|         1.15|
|[kc-c4](kc-c4.c) (2+4 threads) |<=32-mer; <1024 count|            7.5|        35.1|         1.27|
|[yak-count](yak-count.c) (2+4; >=2 count)|<=32-mer; <1024 count| 14.6|        54.8|         0.47|
|[jellyfish2][jf] (16 threads)  |                    |            10.8|       163.9|         0.82|
|[KMC3][KMC] (16 thr; in-mem)   |                    |             9.2|        36.2|         5.02|

## Discussions

* [kc-py1.py](kc-py1.py) is a basic Python3 implementation. It uses string
  translate for fast complementary. Interestingly, pypy is much slower than
  python3. Perhaps the official python3 comes with a better hash table
  implementation. Just a guess. I often recommend pypy over python. I need to
  be more careful about this recommendation in future.

* [kc-cpp1.cpp](kc-cpp1.cpp) implements a basic counter in C++11 using STL's
  [unordered\_map][unordermap]. It is slower than python3. This is partly
  because STL's hash table implementation is very inefficient. C++ does not
  necessarily lead to a fast implementation.

* [kc-cpp2.cpp](kc-cpp2.cpp) replaces `std::unordered_map` with Martin Ankerl's
  [robin\_hood][rhhash] hash table library, which is [among the
  fastest][rhbench] hash table implementations. It is now faster than
  kc-py1.py, though the performance gap is small.

* [kc-c1.c](kc-c1.c) packs k-mers no longer than 32bp into 64-bit integers.
  This dramatically improves speed and reduces the peak memory. Most practical
  k-mer counters employs bit packing. Excluding library files, this counter has
  less than 100 coding lines, not much more complex than the C++ or the python
  implementations.

* [kc-c2.c](kc-c2.c) uses an ensemble of hash tables to save 8 bits for
  counter. This reduces the peak memory. The key advantage of using multiple
  hash tables is to implement multithreading. See below.

* [kc-c3.c](kc-c3.c) puts file reading and parsing into a separate thread. The
  performance improvement is minor here, but it sets the stage for the next
  multi-threaded implementation.

* [kc-c4.c](kc-c4.c) is the fastest counter in this series. Due to the use of
  an ensembl of hash tables in kc-c2, we can parallelize the insertion of a
  batch of k-mers. It is much faster than the previous versions. Notably, kc-c4
  also uses less CPU time. This is probably because batching helps data
  locality.

* [yak-count.c](yak-count.c) is adapted from [yak][yak] and uses the same kc-c4
  algorithm. Similar to [BFCounter][BFCnt], it optionally adds a bloom filter
  to filter out most singleton k-mers (k-mers occurring only once in the
  input). Yak needs to update the bloom filter, read the input twice and count
  twice. It is slower but uses less memory. Yak-count is the most complex
  example in this repo, but it is still short. Its code is also better
  organized. Command line: `-b30` (bloom filter with 1 billion bits).

* [jellyfish2][jf] is probably the fastest in-memory k-mer counter to date. It
  uses less memory and more flexible than kc-c4, but it is slower and much more
  complex. Command line: `count -m 31 -C -s 100000000 -o /dev/null -t 16`.

* [KMC3][KMC] is one of the fastest k-mer counters. It uses minimizers and
  relies on sorting. KMC3 is run in the in-memory mode here. The disk mode is
  as fast. KMC3 is optimized for counting much larger datasets. Although it
  uses more RAM here, it generally uses less RAM than jellyfish2 and other
  in-memory counters given high-coverage human data. Command line: `-k31 -t16
  -r -fa`.

## Conclusions

The k-mer counters here are fairly basic implementations only using generic
hash tables. Nonetheless, we show better engineering can carry the basic idea a
long way. If you want to implement your own k-mer counter,
[yak-count.c](yak-count.c) could be a good starting point. It is fast and
relatively simple. By the way, if you have an efficient and simple k-mer
counter (implemented in a few files), please let me know. I will be happy to add it to the table.

[jf]: http://www.genome.umd.edu/jellyfish.html
[unordermap]: http://www.cplusplus.com/reference/unordered_map/unordered_map/
[rhhash]: https://github.com/martinus/robin-hood-hashing
[rhbench]: https://martin.ankerl.com/2019/04/01/hashmap-benchmarks-01-overview/
[gage-b]: https://ccb.jhu.edu/gage_b/datasets/index.html
[yak]: https://github.com/lh3/yak
[BFCnt]: https://github.com/pmelsted/BFCounter
[KMC]: https://github.com/refresh-bio/KMC
