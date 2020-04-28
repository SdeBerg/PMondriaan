#include <vector>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/contraction.hpp"
#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

void contraction::merge_free_vertices(bulk::world& world, pmondriaan::hypergraph& H) {
    auto total_weight = remove_free_vertices_(H);
    recompute_global_size(world, H);
    local_free_weight_ = total_weight;
    global_free_weight_ = bulk::sum(world, total_weight);
}

void contraction::merge_free_vertices(pmondriaan::hypergraph& H) {
    auto total_weight = remove_free_vertices_(H);
    H.set_global_size(H.size());
    local_free_weight_ = total_weight;
    global_free_weight_ = total_weight;
}

/**
 * Assign  the free vertices greedily to optimize the weight balance. Returns
 * the weights of the parts.
 */
std::vector<long> contraction::assign_free_vertices(pmondriaan::hypergraph& H,
                                                    long max_weight_0,
                                                    long max_weight_1,
                                                    std::mt19937& rng) {
    auto weight_parts = H.weight_all_parts(2);
    assign_free_vertices_(H, weight_parts, max_weight_0, max_weight_1, rng);
    return weight_parts;
}

/**
 * Assign the free vertices of a parallel hypergraph greedily. Returns
 * the weights of the parts.
 */
std::vector<long> contraction::assign_free_vertices(bulk::world& world,
                                                    pmondriaan::hypergraph& H,
                                                    long max_weight_0,
                                                    long max_weight_1,
                                                    std::mt19937& rng) {
    auto weight_parts = global_weight_parts(world, H, 2);

    auto local_free_weights = bulk::coarray<long>(world, world.active_processors());
    for (auto t = 0; t < world.active_processors(); t++) {
        local_free_weights(t)[world.rank()] = local_free_weight_;
    }
    world.sync();
    long total_free_weight = 0;
    for (auto t = 0; t < world.active_processors(); t++) {
        total_free_weight += local_free_weights[t];
    }

    long extra_weight_0 = max_weight_0 - weight_parts[0];
    long extra_weight_1 = max_weight_1 - weight_parts[1];
    long assign_to_part_0 =
    (extra_weight_0 * total_free_weight) / (extra_weight_0 + extra_weight_1);
    // We find the processor whose free vertices have to be split over the parts
    long assigned = 0;
    int split_proc = 0;
    while ((split_proc < world.active_processors()) &&
           (assigned + local_free_weights[split_proc] <= assign_to_part_0)) {
        assigned += local_free_weights[split_proc];
        split_proc++;
    }
    auto new_weight_0 = bulk::var<long>(world);
    auto new_weight_1 = bulk::var<long>(world);

    if (world.rank() < split_proc) {
        assign_all_vertices_(H, 0);
    } else if (world.rank() > split_proc) {
        assign_all_vertices_(H, 1);
    } else {
        // We compute the new part weights assigned so far
        weight_parts[0] += assigned;
        for (auto t = split_proc + 1; t < world.active_processors(); t++) {
            weight_parts[1] += local_free_weights[t];
        }
        assign_free_vertices_(H, weight_parts, max_weight_0, max_weight_1, rng);
        new_weight_0.broadcast(weight_parts[0]);
        new_weight_1.broadcast(weight_parts[1]);
    }
    world.sync();
    weight_parts[0] = new_weight_0;
    weight_parts[1] = new_weight_1;
    return weight_parts;
}

long contraction::remove_free_vertices_(pmondriaan::hypergraph& H) {
    long total_weight = 0;
    for (auto& v : H.vertices()) {
        if (v.degree() == 0) {
            add_free_vertex_(v.id(), v.weight());
            total_weight += v.weight();
        }
    }
    for (auto free_vertex : free_vertices_) {
        H.remove_free_vertex(free_vertex.first);
    }
    return total_weight;
}

void contraction::assign_all_vertices_(pmondriaan::hypergraph& H, long part) {
    for (auto free_vertex : free_vertices_) {
        auto id = free_vertex.first;
        auto weight = free_vertex.second;
        H.add_vertex(id, std::vector<long>(), weight);
        H(id).set_part(part);
    }
}

/**
 * Assign  the free vertices greedily to optimize the weight balance.
 */
void contraction::assign_free_vertices_(pmondriaan::hypergraph& H,
                                        std::vector<long>& weight_parts,
                                        long max_weight_0,
                                        long max_weight_1,
                                        std::mt19937& rng) {
    std::sort(free_vertices_.begin(), free_vertices_.end(),
              [](auto lhs, auto rhs) { return lhs.second > rhs.second; });
    for (auto free_vertex : free_vertices_) {
        auto id = free_vertex.first;
        auto weight = free_vertex.second;
        long part = 0;
        if ((weight_parts[0] + weight > max_weight_0) &&
            (weight_parts[1] + weight > max_weight_1)) {
            std::cout
            << "Could not make a balanced partition using free vertices\n";
        } else if (weight_parts[0] + weight > max_weight_0) {
            part = 1;
        } else if (weight_parts[1] + weight > max_weight_1) {
            part = 0;
        } else if (weight_parts[0] * max_weight_1 < weight_parts[1] * max_weight_0) {
            part = 0;
        } else if (weight_parts[0] * max_weight_1 > weight_parts[1] * max_weight_0) {
            part = 1;
        } else {
            part = rng() % 2;
        }
        weight_parts[part] += weight;
        H.add_vertex(id, std::vector<long>(), weight);
        H(id).set_part(part);
    }
}

} // namespace pmondriaan
