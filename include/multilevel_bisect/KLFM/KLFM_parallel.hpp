#pragma once

#include <limits>
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
              long cut_size = std::numeric_limits<long>::max());

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
 * Initializes the previous_C counts using communication and return the cutsize
 * of the nets the processor is responsible for.
 */
long init_previous_C(bulk::world& world,
                     pmondriaan::hypergraph& H,
                     std::vector<std::vector<long>>& C,
                     bulk::coarray<long>& previous_C,
                     bulk::coarray<long>& cost_my_nets,
                     bulk::block_partitioning<1>& net_partition);

/**
 * Updates the gain values that were outdated.
 */
void update_gains(pmondriaan::hypergraph& H,
                  pmondriaan::net& net,
                  std::vector<long> C_loc,
                  std::vector<long> C_new,
                  pmondriaan::gain_structure& gain_structure);

/**
 * Finds the best moves for a processor sequentially, by only updating local data.
 */
void find_top_moves(pmondriaan::hypergraph& H,
                    pmondriaan::gain_structure& gain_structure,
                    std::vector<std::tuple<long, long, long>>& moves,
                    std::array<long, 2>& weights,
                    long max_weight_0,
                    long max_weight_1,
                    std::mt19937& rng);

/**
 * Determines how many moves from part 0 or 1 should be rejected on each processor.
 * Return the number of rejected moves for each processor. This is positive when
 * we have to move back vertices from part 0 to part 1 and positive otherwise.
 */
long reject_unbalanced_moves(bulk::world& world,
                             bulk::queue<long, long, int, long>& moves_queue,
                             std::array<long, 2>& total_weights,
                             long max_weight_0,
                             long max_weight_1);

/**
 * Updates the counts by communicating with the responsible processor, returns the new cutsize_my_nets.
 */
long update_C(bulk::world& world,
              pmondriaan::hypergraph& H,
              std::vector<std::vector<long>>& C,
              bulk::coarray<long>& previous_C,
              bulk::queue<long, long>& update_nets,
              bulk::partitioning<1>& net_partition,
              bulk::coarray<long>& cost_my_nets,
              long cut_size_my_nets,
              bool update_g,
              pmondriaan::gain_structure& gain_structure);

// For testing purposes
void check_C(bulk::world& world, pmondriaan::hypergraph& H, std::vector<std::vector<long>>& C);

} // namespace pmondriaan
