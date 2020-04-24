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
    global_free_weight_ = bulk::sum(world, total_weight);
}

void contraction::merge_free_vertices(pmondriaan::hypergraph& H) {
    auto total_weight = remove_free_vertices_(H);
    H.set_global_size(H.size());
    global_free_weight_ = total_weight;
}

long contraction::remove_free_vertices_(pmondriaan::hypergraph& H) {
    long total_weight = 0;
    for (auto& v : H.vertices()) {
        if (v.degree() == 0) {
            std::cout << "Add free vertex " << v.id() << "\n";
            add_free_vertex_(v.id(), v.weight());
            total_weight += v.weight();
        }
    }
    for (auto free_vertex : free_vertices_) {
        H.remove_free_vertex(free_vertex.first);
    }
    return total_weight;
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
    return weight_parts;
}

} // namespace pmondriaan
