#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "KLFM.hpp"
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
          int max_passes) {
    int pass = 0;
    bool improvement = true;
    long quality_prev_solution = LONG_MAX;

    while ((pass < max_passes) && improvement) {
        auto result =
        KLFM_pass(world, H, quality_prev_solution, max_weight_0, max_weight_1, labels);
        if (result < quality_prev_solution) {
            quality_prev_solution = result;
        } else {
            improvement = false;
        }
    }
}

/**
 * Runs a single pass of the KLFM algorithm to improve a given partitioning.
 */
long KLFM_pass(bulk::world& world,
               pmondriaan::hypergraph& H,
               long solution_quality,
               long max_weight_0,
               long max_weight_1,
               interval labels) {

    // init gain data structure, also make sure this adhere to the max weights

    long best_solution_quality = solution_quality;
    auto no_improvement_moves = std::vector<int>();

    while (!gain_structure.done()) {
        solution_quality += gain_structure.gain_next();
        auto v_to_move = gain_structure.next();
        gain_structure.update_gains();
        lock(v_to_move); // or is this already okay if we do not enter it into the other gain queue?
        H.move(v_to_move, labels);

        if (solution_quality > best_solution_quality) {
            no_improvement_moves.push_back(v_to_move);
        } else {
            no_improvement_moves.clear();
        }
    }

    // rollback until we are at the best solution seen this pass
    for (auto v : no_improvement_moves) {
        H.move(v, labels);
    }

    return best_solution_quality;
}

} // namespace pmondriaan
