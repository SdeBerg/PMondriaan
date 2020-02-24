#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"
#include <algorithm.hpp>

namespace pmondriaan {

void vertex::remove_net(int n) {
    auto it = std::find(nets_.begin(), nets_.end(), n);
    if (it != nets_.end()) {
        std::iter_swap(it, nets_.end() - 1);
        nets_.erase(nets_.end() - 1);
    }
}

long hypergraph::total_weight() {
    long total = 0;
    for (auto v : vertices_) {
        total += v.weight();
    }
    return total;
}

// computes the sum of the weights of vertices in part
long hypergraph::weight_part(int part) {
    long total = 0;
    for (auto v : vertices_) {
        if (v.part() == part) {
            total += v.weight();
        }
    }
    return total;
}

// computes the weights of all parts upto kbhit
std::vector<long> hypergraph::weight_all_parts(int k) {
    auto total = std::vector<long>(k);
    for (auto v : vertices_) {
        total[v.part()] += v.weight();
    }
    return total;
}

void hypergraph::add_to_nets(pmondriaan::vertex& v) {
    for (auto net_id : v.nets()) {
        nets_[net_id].vertices().push_back(v.id());
    }
}

// removes id from all nets
void hypergraph::remove_from_nets(int id) {
    for (auto& net : nets_) {
        auto it = std::find(net.vertices().begin(), net.vertices().end(), id);
        if (it != net.vertices().end()) {
            std::iter_swap(it, net.vertices().end() - 1);
            net.vertices().erase(net.vertices().end() - 1);
        }
    }
}

void hypergraph::update_map() {
    for (auto i = 0u; i < vertices_.size(); i++) {
        global_to_local[vertices_[i].id()] = i;
    }
}

void hypergraph::set_global_net_sizes(std::vector<size_t>& sizes) {
    for (auto i = 0u; i < nets_.size(); i++) {
        nets_[i].set_global_size(sizes[i]);
    }
}

void hypergraph::print() {

    for (auto& v : vertices_) {
        std::cout << v.id() << ": ";
        for (auto n : v.nets()) {
            std::cout << n << " ";
        }
        std::cout << "\n";
    }
}

/**
 * Compute the global weight of a hypergraph.
 */
long global_weight(bulk::world& world, pmondriaan::hypergraph& H) {

    bulk::var<long> local_weight(world);
    local_weight = H.total_weight();

    long global_weight =
    bulk::foldl(local_weight, [](auto& lhs, auto rhs) { lhs += rhs; });
    return global_weight;
}

/**
 * Compute the global load imbalance of a hypergraph split into k parts.
 */
double load_balance(bulk::world& world, pmondriaan::hypergraph& H, int k) {

    auto weight_parts_coar = bulk::coarray<long>(world, k);
    auto weight_parts = H.weight_all_parts(k);

    for (int i = 0; i < k; i++) {
        weight_parts_coar[i] = weight_parts[i];
    }

    long global_weight = pmondriaan::global_weight(world, H);

    // compute the global part weights
    weight_parts = pmondriaan::foldl_each(weight_parts_coar,
                                          [](auto& lhs, auto rhs) { lhs += rhs; });

    // we compute the global part with largest weight
    long max_weight_part = *std::max_element(weight_parts.begin(), weight_parts.end());

    double eps = ((double)(max_weight_part * k) / (double)global_weight) - 1.0;

    return eps;
}

/**
 * Compute the cutsize with the correct metric
 */
long cutsize(bulk::world& world, pmondriaan::hypergraph& H, int k, std::string metric) {

    long result = 0;
    // this queue contains all labels present for each net
    auto labels = bulk::queue<int, int[]>(world);
    for (auto& net : H.nets()) {
        if (net.size() > 0) {
            auto labels_net = std::set<int>();
            for (auto& v : net.vertices()) {
                labels_net.insert(H(H.local_id(v)).part());
            }
            for (int t = 0; t < world.active_processors(); t++) {
                labels(t).send(net.id(), std::vector<int>(labels_net.begin(),
                                                          labels_net.end()));
            }
        }
    }

    world.sync();

    auto total_cut = std::vector<std::set<int>>(H.nets().size(), std::set<int>());
    for (const auto& [net, labels_net] : labels) {
        for (auto l : labels_net) {
            total_cut[net].insert(l);
        }
    }

    for (auto i = 0u; i < H.nets().size(); i++) {
        if (metric == "cutnet") {
            if (total_cut[i].size() > 1) {
                result += H.net(i).cost();
            }
        } else {
            if (metric == "lambda1") {
                if (total_cut[i].size() > 1) {
                    result += (total_cut[i].size() - 1) * H.net(i).cost();
                }
            } else {
                std::cerr << "Error: unknown metric";
                break;
            }
        }
    }

    return result;
}


/**
 * Compute the global net sizes of a hypergraph.
 */
std::vector<size_t> global_net_sizes(bulk::world& world, pmondriaan::hypergraph& H) {
    auto& nets = H.nets();
    auto net_sizes_coar = bulk::coarray<size_t>(world, nets.size());

    for (auto i = 0u; i < nets.size(); i++) {
        net_sizes_coar[i] = nets[i].size();
    }

    auto net_sizes = pmondriaan::foldl_each(net_sizes_coar, [](auto& lhs, auto rhs) {
        lhs += rhs;
    });

    H.set_global_net_sizes(net_sizes);

    return net_sizes;
}


/**
 * Removes all free nets.
 */
void remove_free_nets(bulk::world& world, pmondriaan::hypergraph& H) {

    auto net_sizes = global_net_sizes(world, H);

    for (auto n = 0u; n < H.nets().size(); n++) {
        if (net_sizes[n] == 1) {
            H.net(n).set_global_size(0);
            if (H.net(n).size() == 1) {
                int v = H.net(n).vertices()[0];
                H(H.local_id(v)).remove_net(n);
                H.net(n).vertices().clear();
            }
        }
    }
}

/**
 * Creates a new hypergraph that only contains the vertices of H with local id between start and end.
 */
pmondriaan::hypergraph
create_new_hypergraph(bulk::world& new_world, pmondriaan::hypergraph& H, int start, int end) {

    std::vector<pmondriaan::vertex> new_vertices(H.vertices().begin() + start,
                                                 H.vertices().begin() + end);
    auto new_nets = std::vector<pmondriaan::net>();
    for (auto& n : H.nets()) {
        new_nets.push_back(pmondriaan::net(n.id(), std::vector<int>()));
    }

    for (auto& v : new_vertices) {
        for (auto n : v.nets()) {
            new_nets[n].add_vertex(v.id());
        }
    }

    bulk::var<int> new_size(new_world);
    new_size = (int)new_vertices.size();
    auto new_global_size =
    bulk::foldl(new_size, [](int& lhs, int rhs) { lhs += rhs; });

    auto new_H = pmondriaan::hypergraph(new_global_size, new_vertices, new_nets);
    remove_free_nets(new_world, new_H);

    return new_H;
}


} // namespace pmondriaan