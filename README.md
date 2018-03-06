Barnes2018-Landscape
============================

**Title of Manuscript**:
Accelerating a landscape evolution model with parallelism

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

Solving inverse problems and achieving statistical rigour in landscape evolution
models requires running many model realizations. Parallel computation is
necessary to achieve this in a reasonable time. However, no previous algorithm
is well-suited to leveraging modern parallelism. Here, I describe an algorithm
that can utilize the parallel potential of GPUs, many-core processors, and SIMD
instructions, in addition to working well in serial. The new algorithm runs 43 x
faster (70 s vs. 3,000 s on a 10,000 x 10,000 input) than the previous state of
the art and exhibits sublinear scaling with input size. I also identify methods
for using multidirectional flow routing and quickly eliminating landscape
depressions and local minima. Tips for parallelization and a step-by-step guide
to achieving it are given to help others achieve good performance with their own
code. Complete, well-commented, easily adaptable source code for all versions of
the algorithm is available as a supplement and on Github.



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
