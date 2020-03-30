#include <limits.h>
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
 * Runs the KLFM algorithm in parallel to improve a given partitioning. Return the quality of the best solution.
 */
long KLFM_par(bulk::world& world,
              pmondriaan::hypergraph& H,
              std::vector<std::vector<long>>& C,
              long weight_0,
              long weight_1,
              long max_weight_0,
              long max_weight_1,
              pmondriaan::options& opts,
              std::mt19937& rng,
              long cut_size = LONG_MAX);

/**
 * Runs a single pass of the KLFM algorithm to improve a given partitioning.
 */
long KLFM_pass_par(bulk::world& world,
                   pmondriaan::hypergraph& H,
                   std::vector<std::vector<long>>& C,
                   long cut_size,
                   std::array<long, 2>& weights,
                   long max_weight_0,
                   long max_weight_1,
                   std::mt19937& rng);

void find_top_moves(pmondriaan::hypergraph& H,
                    pmondriaan::gain_structure& gain_structure,
                    std::vector<std::pair<int, long>> moves,
                    std::array<long, 2>& weights,
                    long max_weight_0,
                    long max_weight_1,
                    std::mt19937& rng);

} // namespace pmondriaan
