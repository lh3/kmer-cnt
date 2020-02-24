## Getting Started

```sh
git clone https://github.com/lh3/kmer-cnt
cd kmer-cnt
make  # C++11 required to compile the two C++ implementations
wget https://github.com/lh3/kmer-cnt/releases/download/v0.1/M_abscessus_HiSeq_10M.fa.gz
./kc-c4 M_abscessus_HiSeq_10M.fa.gz > kc-c4.out
```

## Introduction

K-mer counting is the foundation of many mappers, assemblers and miscellaneous
tools (e.g. genotypers, metagenomics profilers, etc). It is one of the most
important algorithms in Bioinformatics. Here we will implement basic k-mer
counting algorithms but with advanced engineering tricks. We will see how far
better engineering can go.

## Results

We provide seven k-mer counters, which are detailed below the result table. All
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
|[kc-c1](kc-c1.c)               |<=32-mer            |            39.0|        39.0|         1.52|
|[kc-c2](kc-c2.c)               |<=32-mer; <256 count|            37.8|        37.8|         1.02|
|[kc-c3](kc-c3.c)               |<=32-mer; <256 count|            36.0|        40.6|         1.07|
|[kc-c4](kc-c4.c) (2+4 threads) |<=32-mer; <256 count|             8.4|        37.6|         1.24|
|[jellyfish2][jf] (16 threads)  |                    |            13.6|       212.7|         0.82|

Among these k-mer counters:

* [kc-py1.py](kc-py1.py) is a basic Python3 implementation. It uses string
  translate for fast complementary. Interestingly, pypy is much slower than
  python3. Perhaps the official python3 comes with a better hash table
  implementation. Just a guess.

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
  k-mer counters employs bit packing.

* [kc-c2.c](kc-c2.c) uses a collection of hash tables to save 8 bits for
  counter. This reduces the peak memory.

* [kc-c3.c](kc-c3.c) puts file reading and parsing into a separate thread. The
  performance improvement is minor here, but it sets the stage for the next
  multi-threaded implementation.

* [kc-c4.c](kc-c4.c) is the fastest counter in this series. Due to the use of a
  collection of hash tables in kc-c2, we can parallelize the insertion of a
  batch of k-mers. On this particular test dataset, kc-c4 is faster than
  the fastest in-memory k-mer counter [jellyfish2][jf] to date. Although kc-c4
  uses more memory and imposes more restrictions, it is much simpler.

## Conclusions

Better engineering may make big impact on performance.

[jf]: http://www.genome.umd.edu/jellyfish.html
[unordermap]: http://www.cplusplus.com/reference/unordered_map/unordered_map/
[rhhash]: https://github.com/martinus/robin-hood-hashing
[rhbench]: https://martin.ankerl.com/2019/04/01/hashmap-benchmarks-01-overview/
[gage-b]: https://ccb.jhu.edu/gage_b/datasets/index.html
