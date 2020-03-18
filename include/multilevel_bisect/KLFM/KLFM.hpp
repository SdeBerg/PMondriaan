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
          interval labels,
          int max_passes);

/**
 * Runs a single pass of the KLFM algorithm to improve a given partitioning.
 */
long KLFM_pass(bulk::world& world,
               pmondriaan::hypergraph& H,
               long solution_quality,
               long max_weight_0,
               long max_weight_1,
               interval labels);

} // namespace pmondriaan
