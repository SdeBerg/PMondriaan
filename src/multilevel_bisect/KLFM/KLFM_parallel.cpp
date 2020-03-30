#include <array>
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
 * Runs the KLFM algorithm in parallel to improve a given partitioning. Return the cutsize of the best solution.
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
              long cut_size) {

    int pass = 0;
    long prev_cut_size;
    if (cut_size == LONG_MAX) {
        prev_cut_size = pmondriaan::cutsize(world, H, opts.metric);
    } else {
        prev_cut_size = cut_size;
    }
    auto total_weights = std::array<long, 2>();
    total_weights[0] = weight_0;
    total_weights[1] = weight_1;
    while (pass < opts.KLFM_max_passes) {
        auto result = KLFM_pass_par(world, H, C, prev_cut_size, total_weights,
                                    max_weight_0, max_weight_1, rng);
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
 * Runs a single pass of the KLFM algorithm in parallel to improve a given partitioning.
 */
long KLFM_pass_par(bulk::world& world,
                   pmondriaan::hypergraph& H,
                   std::vector<std::vector<long>>& C,
                   long cut_size,
                   std::array<long, 2>& total_weights,
                   long max_weight_0,
                   long max_weight_1,
                   std::mt19937& rng) {

    auto max_extra_weight = std::array<long, 2>();
    auto gain_structure = pmondriaan::gain_structure(H, C);

    long best_cut_size = cut_size;
    auto no_improvement_moves = std::vector<int>();

    while (not everyone done) {
        // find best KLFM_par_number_send_moves moves
        auto moves = std::vector < std::pair<int, long>(opts.KLFM_par_number_send_moves);
        find_top_moves(H, gain_structure, moves, total_weights, max_weight_0,
                       max_weight_1, rng);

      send moves to processor 0
      if (world.rank() == 0) {
        select the top moves that adhere to the balance constraint
        inform procs how many have been selected and where no gain moves start
        send new balance
        send whether or not to reset no gain moves
      }
      update balance
      reverse moves that have not been selected
      send updates about counts in nets to responsible processors
      send updated counts of nets you are responsible for to processors that need it
      update C and update gains
    }
    reverse moves tot beste punt

    /*while (!gain_structure.done()) {
        max_extra_weight[0] = max_weight_0 - weights[0];
        max_extra_weight[1] = max_weight_1 - weights[1];

        int part_to_move =
        gain_structure.part_next(max_extra_weight[0], max_extra_weight[1], rng);
        auto v_to_move = gain_structure.next(part_to_move);

        if (max_extra_weight[(part_to_move + 1) % 2] - H(H.local_id(v_to_move)).weight() >= 0) {

            cut_size -= gain_structure.gain_next(part_to_move);
            weights[part_to_move] -= H(H.local_id(v_to_move)).weight();
            weights[(part_to_move + 1) % 2] += H(H.local_id(v_to_move)).weight();
            // std::cout << "weights: " << weights[0] << " " << weights[1] << "\n";
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
    }*/

    return best_cut_size;
}

void find_top_moves(pmondriaan::hypergraph& H,
                    pmondriaan::gain_structure& gain_structure,
                    std::vector<std::tuple<int, long, long>> moves,
                    std::array<long, 2>& weights,
                    long max_weight_0,
                    long max_weight_1,
                    std::mt19937& rng) {
    int moves_found = 0;
    while (!gain_structure.done() && (moves_found < moves)) {
        max_extra_weight[0] = max_weight_0 - weights[0];
        max_extra_weight[1] = max_weight_1 - weights[1];

        int part_to_move =
        gain_structure.part_next(max_extra_weight[0], max_extra_weight[1], rng);
        auto v_to_move = gain_structure.next(part_to_move);

        if (max_extra_weight[(part_to_move + 1) % 2] - H(H.local_id(v_to_move)).weight() >= 0) {
            auto weight_v = H(H.local_id(v_to_move)).weight();
            if (part_to_move == 0) {
                moves[moves_found] =
                std::make_tuple(v_to_move, gain_structure.gain_next(part_to_move),
                                -1 * weight_v);
                weights[0] -= weight_v;
                weights[1] += weight_v;
            } else {
                moves[moves_found] =
                std::make_tuple(v_to_move, gain_structure.gain_next(part_to_move), weight_v);
                weights[1] -= weight_v;
                weights[0] += weight_v;
            }
            gain_structure.move(v_to_move);
        } else {
            gain_structure.remove(v_to_move);
        }
    }
}

} // namespace pmondriaan
