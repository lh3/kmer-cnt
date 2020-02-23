# The Making of a Faster k-mer Counter

This README is organized as a blog post. The impatient may jump to the end of
the README to find the conclusions.

## Conclusion

|Implementation                 | Elapsed time (s) | CPU time (s) | Peak RAM (GB) |
|:------------------------------|-----------------:|-------------:|--------------:|
|[kc-py1](kc-py1.py) + Python3.7|             499.6|         499.5|           8.15|
|[kc-py1](kc-py1.py) + Pypy7.3  |            1220.8|        1220.8|          12.21|
|[kc-cpp1](kc-cpp1.cpp)         |             528.0|         527.9|           8.27| 
|[kc-cpp2](kc-cpp2.cpp)         |             319.6|         319.6|           6.90|
|[kc-c1](kc-c1.c) (<=32-mer)       |           39.0|          39.0|           1.52|
|[kc-c2](kc-c2.c) (<=32-mer; <=255)|           37.8|          37.8|           1.02|
|[kc-c3](kc-c3.c) (<=32-mer; <=255)|           36.0|          40.6|           1.07|
|[kc-c4](kc-c4.c) (<=32-mer; <=255)|            8.4|          37.6|           1.24|
|[jellyfish2][jf] (16 threads)     |           13.6|         212.7|           0.82|

[jf]: http://www.genome.umd.edu/jellyfish.html
