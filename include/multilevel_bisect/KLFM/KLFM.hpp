#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

/**
 * Runs the KLFM algorithm to improve a given partitioning.
 */
void KLFM(bulk::world& world,
                          pmondriaan::hypergraph& H,
                          long max_weight_0,
                          long max_weight_1,
                          interval labels);

} // namespace pmondriaan
