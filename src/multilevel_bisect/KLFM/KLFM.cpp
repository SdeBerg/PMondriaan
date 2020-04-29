#include <array>
#include <limits>
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
 * Runs the KLFM algorithm to improve a given partitioning. Return the cutsize of the best solution.
 */
long KLFM(pmondriaan::hypergraph& H,
          std::vector<std::vector<long>>& C,
          long weight_0,
          long weight_1,
          long max_weight_0,
          long max_weight_1,
          pmondriaan::options& opts,
          std::mt19937& rng,
          long cut_size) {

    size_t pass = 0;
    long prev_cut_size;
    if (cut_size == std::numeric_limits<long>::max()) {
        prev_cut_size = pmondriaan::cutsize(H, opts.metric);
    } else {
        prev_cut_size = cut_size;
    }
    auto weights = std::array<long, 2>();
    weights[0] = weight_0;
    weights[1] = weight_1;
    while (pass < opts.KLFM_max_passes) {
        auto result =
        KLFM_pass(H, C, prev_cut_size, weights, max_weight_0, max_weight_1, rng);
        if (result < prev_cut_size) {
            prev_cut_size = result;
        } else {
            break;
        }
        pass++;
    }
    return prev_cut_size;
}

/**
 * Runs a single pass of the KLFM algorithm to improve a given partitioning.
 */
long KLFM_pass(pmondriaan::hypergraph& H,
               std::vector<std::vector<long>>& C,
               long cut_size,
               std::array<long, 2>& weights,
               long max_weight_0,
               long max_weight_1,
               std::mt19937& rng) {

    auto max_extra_weight = std::array<long, 2>();
    auto gain_structure = pmondriaan::gain_structure(H, C);

    long best_cut_size = cut_size;
    auto no_improvement_moves = std::vector<long>();

    while (!gain_structure.done()) {
        max_extra_weight[0] = max_weight_0 - weights[0];
        max_extra_weight[1] = max_weight_1 - weights[1];

        long part_to_move =
        gain_structure.part_next(max_extra_weight[0], max_extra_weight[1], rng);
        auto v_to_move = gain_structure.next(part_to_move);

        if (max_extra_weight[(part_to_move + 1) % 2] - H(H.local_id(v_to_move)).weight() >= 0) {

            cut_size -= gain_structure.gain_next(part_to_move);
            weights[part_to_move] -= H(H.local_id(v_to_move)).weight();
            weights[(part_to_move + 1) % 2] += H(H.local_id(v_to_move)).weight();
            gain_structure.move(v_to_move);
            if (cut_size > best_cut_size) {
                no_improvement_moves.push_back(v_to_move);
            } else {
                best_cut_size = cut_size;
                no_improvement_moves.clear();
            }
        } else {
            gain_structure.remove(v_to_move);
        }
    }

    // rollback until we are at the best solution seen this pass
    for (auto v : no_improvement_moves) {
        auto& vertex = H(H.local_id(v));
        weights[vertex.part()] -= vertex.weight();
        weights[(vertex.part() + 1) % 2] += vertex.weight();
        H.move(v, C);
    }

    return best_cut_size;
}

} // namespace pmondriaan
