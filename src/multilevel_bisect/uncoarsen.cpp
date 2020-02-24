#include <vector>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/contraction.hpp"
#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/uncoarsen.hpp"

namespace pmondriaan {

/**
 * Uncoarsens the hypergraph HC into the hypergraph H.
 */
void uncoarsen_hypergraph(bulk::world& world,
                          pmondriaan::hypergraph& HC,
                          pmondriaan::hypergraph& H,
                          pmondriaan::contraction& C) {

    /* The part_queue will contain the vertex and the part is has been assigned to in the uncoarsened hypergraph */
    auto part_queue = bulk::queue<int, int>(world);

    for (auto& v : HC.vertices()) {
        H(H.local_id(v.id())).set_part(v.part());
    }

    for (auto i = 0u; i < C.size(); i++) {
        auto sample = C.id_sample(i);
        auto part = HC(HC.local_id(sample)).part();
        for (auto match : C.matches(i)) {
            part_queue(match.proc()).send(match.id(), part);
        }
    }

    world.sync();

    for (const auto& [id, part] : part_queue) {
        H(H.local_id(id)).set_part(part);
    }
}

} // namespace pmondriaan
