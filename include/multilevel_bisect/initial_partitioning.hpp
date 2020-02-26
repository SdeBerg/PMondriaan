#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"
#include "recursive_bisection.hpp"

namespace pmondriaan {

/**
 * Creates an initial partitioning for hypergraph H.
 */
void initial_partitioning(bulk::world& world,
                          pmondriaan::hypergraph& H,
                          long max_weight_0,
                          long max_weight_1,
                          interval labels);

} // namespace pmondriaan
