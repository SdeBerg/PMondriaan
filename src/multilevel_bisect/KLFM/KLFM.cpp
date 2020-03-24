#include <limits.h>
#include <random>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/KLFM/KLFM.hpp"
#include "multilevel_bisect/KLFM/gain_buckets.hpp"

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
          std::mt19937& rng) {

    int pass = 0;
    long quality_prev_solution = pmondriaan::cutsize(H, opts.metric);

    while (pass < opts.KLFM_max_passes) {
        auto result = KLFM_pass(H, C, quality_prev_solution, weight_0, weight_1,
                                max_weight_0, max_weight_1, rng);
        if (result < quality_prev_solution) {
            quality_prev_solution = result;
        } else {
            break;
        }
        pass++;
    }
    return quality_prev_solution;
}

/**
 * Runs a single pass of the KLFM algorithm to improve a given partitioning.
 */
long KLFM_pass(pmondriaan::hypergraph& H,
               std::vector<std::vector<long>>& C,
               long solution_quality,
               long weight_0,
               long weight_1,
               long max_weight_0,
               long max_weight_1,
               std::mt19937& rng) {

    long max_extra_weight_0 = max_weight_0 - weight_0;
    long max_extra_weight_1 = max_weight_1 - weight_1;

    auto gain_structure = pmondriaan::gain_structure(H, C);

    long best_solution_quality = solution_quality;
    auto no_improvement_moves = std::vector<int>();

    while (!gain_structure.done()) {
        int part_to_move =
        gain_structure.part_next(max_extra_weight_0, max_extra_weight_1, rng);
        std::cout << "gain: " << gain_structure.gain_next(part_to_move) << " v "
                  << gain_structure.next(part_to_move) << "\n";
        solution_quality -= gain_structure.gain_next(part_to_move);
        auto v_to_move = gain_structure.next(part_to_move);
        gain_structure.move(v_to_move);

        if (solution_quality > best_solution_quality) {
            no_improvement_moves.push_back(v_to_move);
        } else {
            best_solution_quality = solution_quality;
            no_improvement_moves.clear();
        }
    }

    // rollback until we are at the best solution seen this pass
    for (auto v : no_improvement_moves) {
        H.move(v, C);
    }

    return best_solution_quality;
}

} // namespace pmondriaan
