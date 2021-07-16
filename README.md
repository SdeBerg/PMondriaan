# PMondriaan

PMondriaan is a parallel hypergraph partitioner that can partition a hypergraph into k parts using p processors. Hypergraph partitioning is the engine of the [Mondriaan](https://webspace.science.uu.nl/~bisse101/Mondriaan/Docs/USERS_GUIDE.html) sparse matrix partitioner. Therefore, a parallel 2D sparse matrix partitioner can be built that uses the PMondriaan hypergraph partitioner, for example by using the medium-grain method, hence the name PMondriaan.

PMondriaan was created in the C++ programming language, using the Bulk software library. Bulk is an interface for writing BSP programs in the C++ language. Using Bulk, we were able to write PMondriaan for distributed-memory, shared-memory, and hybrid systems simultaneously. 

## Getting started

The easiest way to get started is to download the source code from this Github page. PMondriaan uses three submodules: Bulk, CLI11, and googletest. To clone the directory including all submodules you can use `git clone --recurse-submodules`. Otherwise, you can initiate the submodules separately using `git submodule update --init --recursive` from your PMondriaan directory.

PMondriaan requires an up-to-date compiler that supports C++17, e.g. GCC >= 7.0, or Clang >= 4.0. To build PMondriaan using CMake do:

    mkdir build
    cd build
    cmake ..
    make

### How to run PMondriaan

PMondriaan can run on distributed-memory, shared-memory, and hybrid systems. Building the PMondriaan library creates both a thread version, for shared-memory systems, and an MPI version, for distributed-memory and hybrid systems, of the PMondriaan program as `Run_PMondriaan_thread` and `Run_PMondriaan_mpi`. The MPIversion is only created if MPI is available on the system. You can now partition, for example, the Dolphins matrix into two parts using the thread backend on two processors by:

    ./Run_PMondriaan_thread -f ../test/data/matrices/dolphins/dolphins.mtx -p 2 -k 2

### Unit tests
### Examples

## Output

The distributed hypergraph graph is stored using an adapted Matrix Market format, similar to the one used in [Mondriaan](https://webspace.science.uu.nl/~bisse101/Mondriaan/Docs/USERS_GUIDE.html) (using the row-net model to translate the hypergraph into a sparse matrix). The structure of the format is as follows:

```
%%MatrixMarket distributed-matrix coordinate pattern general
m n nnz k
Kstart[0] (this should be 0)
…
…
…
Kstart[k] (this should be nnz)
A.i[0] A.j[0]
…
…
…
A.i[nnz-1] A.j[nnz-1]
```
Here, `Kstart` contains the start indices of the parts. For example, the first `Kstart[1]`nonzeros are assigned to part 0.  Contrary to the format used in Mondriaan, we do not include the values of the nonzeros. These values are not used during computation and are therefore not stored when reading the original matrix.  The `MondriaanStats` tool from the [Mondriaan package](https://webspace.science.uu.nl/~bisse101/Mondriaan/Docs/USERS_GUIDE.html) can be used for the analysis of the partitioning. The cut size and load balance of the partitioning are printed at the end of the run.

## Program options

A number of options and parameters are available to run PMondriaan. These can be set in the `tools/defaults.toml` file, or passed as command line arguments in the program call. Input given via the command line overrules the settings in the defaults file. The default values of the options are indicated below by a \*.

### Nonnumerical options

The nonnumerical options are used to define the hypergraph and to choose partitioning methods.

Option | Values | Description
------------ | ------------- | ------------- 
`file` | `path/to/file.mtx` | Specifies the matrix (in the Matrix Market file format) that should be used to create the hypergraph using the row-net model.
`weights` | `one`, `degree*` | How to set the weights of the vertices. The `one` option sets all vertex weights to 1, the `degree` option sets the weight of each vertex to its degree.
`metric` | `cutnet`, `lambda_minus_one*` | Cut metric to be minimized, either the hyperedge-cut or the lambda-minus-one-cut metric.
`bisect` | `random`, `multilevel*` | Bisection method to be used. The random option is only meant for debugging.
`sampling` | `random*`, `label_propagation` | Sampling method to be used. In the label propagation method, a label propagation step is included to select samples that differ significantly. The random method selects samples uniformly at random.

### Numerical options

The numerical options are often used to optimise given partitioning methods. Some should be adjusted to the problem at hand.

Parameter | Default | Description
------------ | ------------- | ------------- 
`k` | 2 | Integer. Range >= 2. Number of parts to partition the hypergraph into.
`p` | 4 | Integer. Range >= 1. Number of processors to be used.
`eps` | 0.03 | Double. Range = [0,1]. Balance constraint the that the output partitioning into k parts should adhere to.
`eta` | 0.03 | Double. Range = [0,1]. Balance constraint of the division of the hypergraph over the p processors that is used during the computation.
`sample_size` | 10000 | Integer. Range >= 1. Global sample size. Should be adjusted to the problem at hand. Usually 1-2% of the vertices should be selected as a sample.
`max_cluster_size` | 50 | Integer. Range >= 2. Recommended range: 20-50. The maximum size of a cluster in the coarsening phase. Should be adjusted for the problem at hand.
`lp_max_iter` | 25 | Integer. Range >= 1. Maximum number of iterations in label propagation step used in the initial partitioning (and in sampling if the label propagation mode is selected).
`coarsening_nrvertices` | 200 | Integer. Range >= 1. Recommended range: 100-500. Determines when to stop coarsening, as the current number of vertices is small enough.
`coarsening_max_rounds` | 128 | Integer. Range >= 1. The maximum number of coarsenings that may be performed.
`KLFM_max_passes` | 25 | Integer. Range >= 1. Maximum number of passes in FM refinement step.
`KLFM_max_no_gain_moves` | 200 | Integer. Range >= 0. Maximum number of successive no-gain moves in the sequential FM refinement.
`KLFM_par_send_moves` | 20 | Integer. Range >= 1. Number of moves generated by each processor before synchronization in the parallel FM refinement algorithm.
