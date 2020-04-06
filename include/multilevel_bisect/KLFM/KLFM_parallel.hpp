#pragma once

#include <limits.h>
#include <random>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/KLFM/gain_buckets.hpp"

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
                   std::array<long, 2>& total_weights,
                   long max_weight_0,
                   long max_weight_1,
                   pmondriaan::options& opts,
                   std::mt19937& rng);

/**
 * Finds the best moves for a processor sequentially, by only updating local data.
 */
void find_top_moves(pmondriaan::hypergraph& H,
                    pmondriaan::gain_structure& gain_structure,
                    std::vector<std::tuple<int, long, long>>& moves,
                    std::array<long, 2>& weights,
                    long max_weight_0,
                    long max_weight_1,
                    std::mt19937& rng);

/**
 * Determines how many moves from part 0 or 1 should be rejected on each processor.
 * Return the number of rejected moves for each processor. This is positive when
 * we have to move back vertices from part 0 to part 1 and positive otherwise.
 */
std::vector<int> reject_unbalanced_moves(int p,
                                         bulk::queue<long, long, int>& moves_queue,
                                         std::array<long, 2>& total_weights,
                                         long max_weight_0,
                                         long max_weight_1);

} // namespace pmondriaan