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
                          interval labels) {

    //bisect_random(world, H, max_weight_0, max_weight_1, 0, H.size(), labels);
    auto L = label_propagation_bisect(H, 100, max_weight_0, max_weight_1);
    for (auto i = 0u; i < H.size(); i++) {
    	H(i).set_part(L[i]);
    }
}


} // namespace pmondriaan
