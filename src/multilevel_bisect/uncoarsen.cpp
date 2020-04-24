#include <random>
#include <unordered_map>
#include <vector>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/contraction.hpp"
#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/KLFM/KLFM.hpp"
#include "multilevel_bisect/KLFM/KLFM_parallel.hpp"
#include "multilevel_bisect/uncoarsen.hpp"

namespace pmondriaan {

/**
 * Uncoarsens the hypergraph HC into the hypergraph H.
 */
void uncoarsen_hypergraph(bulk::world& world,
                          pmondriaan::hypergraph& HC,
                          pmondriaan::hypergraph& H,
                          pmondriaan::contraction& C) {

    // The part_queue will contain the vertex and the part is has been assigned to in the uncoarsened hypergraph
    auto part_queue = bulk::queue<long, long>(world);

    for (auto& v : HC.vertices()) {
        auto id_map = H.map().find(v.id());
        if (id_map != H.map().end()) {
            H(id_map->second).set_part(v.part());
        }
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

/**
 * Uncoarsens the hypergraph HC sequentially into the hypergraph H.
 * The cutsize is then optimized using the KLFM algorithm. Returns
 * the cutsize of the partitioning found.
 */
long uncoarsen_hypergraph_seq(bulk::world& world,
                              pmondriaan::hypergraph& HC,
                              pmondriaan::hypergraph& H,
                              pmondriaan::contraction& C,
                              pmondriaan::options& opts,
                              long max_weight_0,
                              long max_weight_1,
                              long cut_size,
                              std::mt19937& rng) {
    // We first assign the free vertices of HC greedily such that the imbalance is minimized
    auto new_weights = C.assign_free_vertices(HC, max_weight_0, max_weight_1, rng);
    uncoarsen_hypergraph(world, HC, H, C);
    auto counts = pmondriaan::init_counts(H);
    return KLFM(H, counts, new_weights[0], new_weights[1], max_weight_0,
                max_weight_1, opts, rng);
}

/**
 * Uncoarsens the hypergraph HC into the hypergraph H.
 * The cutsize is then optimized using the parallel KLFM algorithm. Returns
 * the cutsize of the partitioning found.
 */
long uncoarsen_hypergraph_par(bulk::world& world,
                              pmondriaan::hypergraph& HC,
                              pmondriaan::hypergraph& H,
                              pmondriaan::contraction& C,
                              pmondriaan::options& opts,
                              long max_weight_0,
                              long max_weight_1,
                              long cut_size,
                              std::mt19937& rng) {
    // We first assign the free vertices of HC greedily such that the imbalance is minimized
    auto new_weights = C.assign_free_vertices(world, HC, max_weight_0, max_weight_1, rng);
    uncoarsen_hypergraph(world, HC, H, C);
    auto counts = pmondriaan::init_counts(world, H);
    return KLFM_par(world, H, counts, new_weights[0], new_weights[1],
                    max_weight_0, max_weight_1, opts, rng, cut_size);
}

} // namespace pmondriaan
