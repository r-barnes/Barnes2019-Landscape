Barnes2018-Landscape
============================

**Title of Manuscript**:
Accelerating a fluvial incision and landscape evolution model

**Authors**: Richard Barnes

**Corresponding Author**: Richard Barnes (richard.barnes@berkeley.edu)

**DOI Number of Manuscript**: TODO

**Code Repositories**
 * [Author's GitHub Repository](https://github.com/r-barnes/Barnes2018-Landscape)

This repository contains a reference implementation of the algorithms presented
in the manuscript above, along with information on acquiring the various
datasets used, and code to perform correctness tests.



Abstract
--------

Solving inverse problems in landscape evolution requires running _many_
model realizations quickly. Braun and Willett (2013) present a collection of
_O(N)_ algorithms for doing so. While their algorithm scales linearly with the
number of cells in a discretized landscape, it is poorly suited to modern
parallelism, contrary to their claims. Here, I describe modifications to the
algorithm that enable it to leverage SIMD for accelerated serial execution;
these same modifications make it suitable for acceleration on GPUs and many-core
processors. The new algorithm runs 43x faster than the original and exhibits
sublinear scaling with input size. I also identify methods for quickly
eliminating landscape depressions and local minima. Complete, well-commented
source code for all versions of the algorithm is available as a supplement and
on Github.



Compilation
-----------

Run

    make

to compile all non-GPU code using your default compiler.

Run

    make -f Makefile.summitdev

to compile all GPU code using the PGI compiler.

Run

    ./tests/tests.sh

To run all tests.
