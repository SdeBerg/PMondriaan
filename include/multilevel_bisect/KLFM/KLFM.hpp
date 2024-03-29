#include <limits>
#include <random>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

/**
 * Runs the KLFM algorithm to improve a given partitioning. Return the quality of the best solution.
 */
long KLFM(pmondriaan::hypergraph& H,
          std::vector<std::vector<long>>& C,
          long weight_0,
          long weight_1,
          long max_weight_0,
          long max_weight_1,
          pmondriaan::options& opts,
          std::mt19937& rng,
          long cut_size = std::numeric_limits<long>::max());

/**
 * Runs a single pass of the KLFM algorithm to improve a given partitioning.
 */
long KLFM_pass(pmondriaan::hypergraph& H,
               std::vector<std::vector<long>>& C,
               long cut_size,
               std::array<long, 2>& weights,
               long max_weight_0,
               long max_weight_1,
               pmondriaan::options& opts,
               std::mt19937& rng);

long make_balanced(pmondriaan::hypergraph& H,
                   std::vector<std::vector<long>>& C,
                   long cut_size,
                   std::array<long, 2>& weights,
                   long max_weight_0,
                   long max_weight_1);

// For testing purposes
void check_C(pmondriaan::hypergraph& H, std::vector<std::vector<long>>& C);

} // namespace pmondriaan
