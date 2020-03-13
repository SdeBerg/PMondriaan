#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "bisect.hpp"
#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/initial_partitioning.hpp"
#include "multilevel_bisect/label_propagation.hpp"

namespace pmondriaan {

/**
 * Creates an initial partitioning for hypergraph H.
 */
void initial_partitioning(bulk::world& world,
                          pmondriaan::hypergraph& H,
                          long max_weight_0,
                          long max_weight_1,
                          interval labels,
                          std::mt19937& rng) {

    // bisect_random(world, H, max_weight_0, max_weight_1, 0, H.size(), labels, rng);
    auto L = label_propagation_bisect(H, 100, max_weight_0, max_weight_1, rng);
    for (auto i = 0u; i < H.size(); i++) {
        if (L[i] == 0) {
            H(i).set_part(labels.low);
        } else {
            H(i).set_part(labels.high);
        }
    }
}


} // namespace pmondriaan
