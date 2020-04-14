#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <vector>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"

#include "algorithm.hpp"
#include "options.hpp"
#include "util/interval.hpp"

namespace pmondriaan {

void vertex::remove_net(int n) {
    auto it = std::find(nets_.begin(), nets_.end(), n);
    if (it != nets_.end()) {
        std::iter_swap(it, nets_.end() - 1);
        nets_.pop_back();
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

// computes the weights of all parts upto k
std::vector<long> hypergraph::weight_all_parts(int k) {
    auto total = std::vector<long>(k);
    for (auto v : vertices_) {
        total[v.part()] += v.weight();
    }
    return total;
}

void hypergraph::add_to_nets(pmondriaan::vertex& v) {
    for (auto net_id : v.nets()) {
        nets_[local_id_net(net_id)].vertices().push_back(v.id());
    }
}

// removes id from all nets
void hypergraph::remove_from_nets(int id) {
    for (auto& net : nets_) {
        auto it = std::find(net.vertices().begin(), net.vertices().end(), id);
        if (it != net.vertices().end()) {
            std::iter_swap(it, net.vertices().end() - 1);
            net.vertices().pop_back();
        }
    }
}

// removes a net and the net from all net lists of vertices
void hypergraph::remove_net_by_index(int index) {
    int id = nets_[index].id();
    for (auto& v : nets_[index].vertices()) {
        vertices_[local_id(v)].remove_net(id);
    }
    std::iter_swap(nets_.begin() + index, nets_.end() - 1);
    nets_.pop_back();
    net_global_to_local.insert_or_assign(nets_[index].id(), index);
    net_global_to_local.erase(id);
}

// moves a vertex to the other part in 0,1
void hypergraph::move(int id) {
    int idl = this->local_id(id);
    auto& vertex = vertices_[idl];
    vertex.set_part((vertex.part() + 1) % 2);
}

void hypergraph::move(int id, std::vector<std::vector<long>>& C) {
    auto& v = vertices_[this->local_id(id)];
    int from = v.part();
    move(id);
    int to = v.part();
    for (auto n : v.nets()) {
        C[local_id_net(n)][from]--;
        C[local_id_net(n)][to]++;
    }
}

void hypergraph::update_map() {
    for (auto i = 0u; i < vertices_.size(); i++) {
        global_to_local[vertices_[i].id()] = i;
    }
}

void hypergraph::update_map_nets() {
    for (auto i = 0u; i < nets_.size(); i++) {
        net_global_to_local[nets_[i].id()] = i;
    }
}

void hypergraph::set_global_net_sizes(std::vector<size_t>& sizes) {
    std::cout << "nets size: " << nets_.size() << " sizes size: " << sizes.size() << "\n";
    for (auto i = 0u; i < nets_.size(); i++) {
        std::cout << "i: " << i << " global size: " << sizes[i] << "\n";
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
 * Initialize the counts for parts 0,1.
 */
std::vector<std::vector<long>> init_counts(pmondriaan::hypergraph& H) {
    auto counts =
    std::vector<std::vector<long>>(H.nets().size(), std::vector<long>(2, 0));
    for (auto& v : H.vertices()) {
        for (auto n : v.nets()) {
            counts[H.local_id_net(n)][v.part()]++;
        }
    }
    return counts;
}

/**
 * Compute the global weight of a hypergraph.
 */
long global_weight(bulk::world& world, pmondriaan::hypergraph& H) {
    auto local_weight = H.total_weight();
    long global_weight = bulk::sum(world, local_weight);
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
    weight_parts =
    bulk::foldl_each(weight_parts_coar, [](auto& lhs, auto rhs) { lhs += rhs; });

    // we compute the global part with largest weight
    long max_weight_part = *std::max_element(weight_parts.begin(), weight_parts.end());

    double eps = ((double)(max_weight_part * k) / (double)global_weight) - 1.0;

    return eps;
}

/**
 * Compute the cutsize with the correct metric of a local hypergraph
 */
long cutsize(pmondriaan::hypergraph& H, pmondriaan::m metric) {
    long result = 0;
    switch (metric) {
    case pmondriaan::m::cut_net: {
        for (auto& net : H.nets()) {
            auto labels_net = std::unordered_set<int>();
            for (auto& v : net.vertices()) {
                labels_net.insert(H(H.local_id(v)).part());
            }
            if (labels_net.size() > 1) {
                result += net.cost();
            }
        }
        break;
    }
    case pmondriaan::m::lambda_minus_one: {
        for (auto& net : H.nets()) {
            auto labels_net = std::unordered_set<int>();
            for (auto& v : net.vertices()) {
                labels_net.insert(H(H.local_id(v)).part());
            }
            if (labels_net.size() > 1) {
                result += (labels_net.size() - 1) * net.cost();
            }
        }
        break;
    }
    default: {
        std::cerr << "Error: unknown metric\n";
        break;
    }
    }
    return result;
}

/**
 * Compute the cutsize with the correct metric
 */
long cutsize(bulk::world& world, pmondriaan::hypergraph& H, pmondriaan::m metric) {
    long result = 0;
    auto net_partition = bulk::block_partitioning<1>({H.global_number_nets()},
                                                     {world.active_processors()});

    // this queue contains all labels present for each net
    auto labels = bulk::queue<int, int[]>(world);
    for (auto& net : H.nets()) {
        auto labels_net = std::unordered_set<int>();
        for (auto& v : net.vertices()) {
            labels_net.insert(H(H.local_id(v)).part());
        }
        labels(net_partition.owner(net.id()))
        .send(net.id(), std::vector<int>(labels_net.begin(), labels_net.end()));
    }
    world.sync();

    auto total_cut =
    std::vector<std::unordered_set<int>>(net_partition.local_size(world.rank())[0]);

    for (const auto& [net, labels_net] : labels) {
        for (auto l : labels_net) {
            total_cut[net_partition.local({net})[0]].insert(l);
        }
    }

    switch (metric) {
    case pmondriaan::m::cut_net: {
        for (auto i = 0u; i < total_cut.size(); i++) {
            if (total_cut[i].size() > 1) {
                result += H.net(i).cost();
            }
        }
        break;
    }
    case pmondriaan::m::lambda_minus_one: {
        for (auto i = 0u; i < total_cut.size(); i++) {
            if (total_cut[i].size() > 1) {
                result += (total_cut[i].size() - 1) * H.net(i).cost();
            }
        }
        break;
    }
    default: {
        std::cerr << "Error: unknown metric\n";
    }
    }
    return bulk::sum(world, result);
}


/**
 * Compute the global net sizes of a hypergraph.
 */
std::vector<size_t> global_net_sizes(bulk::world& world, pmondriaan::hypergraph& H) {

    auto net_partition = bulk::block_partitioning<1>({H.global_number_nets()},
                                                     {world.active_processors()});

    auto nets = H.nets();
    bulk::queue<int, size_t> net_size_queue(world);

    for (auto i = 0u; i < nets.size(); i++) {
        net_size_queue(net_partition.owner(nets[i].id()))
        .send(nets[i].id(), nets[i].size());
    }
    world.sync();

    auto size_nets =
    bulk::coarray<size_t>(world, net_partition.local_count(world.rank()));
    for (auto i = 0; i < net_partition.local_count(world.rank()); i++) {
        size_nets[i] = 0;
    }
    for (const auto& [id, size] : net_size_queue) {
        size_nets[net_partition.local({id})[0]] += size;
    }

    auto result = std::vector<bulk::future<size_t>>(nets.size());
    for (auto i = 0u; i < nets.size(); i++) {
        result[i] =
        size_nets(net_partition.owner(nets[i].id()))[net_partition.local(nets[i].id())[0]]
        .get();
    }
    world.sync();
    std::cout << "net0:" << result[0].value() << "\n ";
    auto r = std::vector<size_t>();
    H.set_global_net_sizes(r);
    return r;
}


/**
 * Removes all free nets.
 */
void remove_free_nets(bulk::world& world, pmondriaan::hypergraph& H) {
    auto net_sizes = global_net_sizes(world, H);
    for (auto n = 0u; n < H.nets().size(); n++) {
        if (net_sizes[n] == 1) {
            H.remove_net_by_index(n);
        }
    }
}

/**
 * Removes all free nets.
 */
void remove_free_nets(pmondriaan::hypergraph& H) {
    for (auto n = 0u; n < H.nets().size(); n++) {
        if (H.net(n).size() <= 1) {
            H.remove_net_by_index(n);
        }
    }
}

/**
 * Creates a new hypergraph that only contains the vertices of H with local id between start and end.
 */
pmondriaan::hypergraph
create_new_hypergraph(bulk::world& new_world, pmondriaan::hypergraph& H, int start, int end) {

    std::cout << "1";
    std::vector<pmondriaan::vertex> new_vertices(H.vertices().begin() + start,
                                                 H.vertices().begin() + end);
    auto new_nets = std::vector<pmondriaan::net>();
    for (auto& n : H.nets()) {
        new_nets.push_back(pmondriaan::net(n.id(), std::vector<int>()));
    }
    std::cout << "2";
    for (auto& v : new_vertices) {
        for (auto n : v.nets()) {
            new_nets[H.local_id_net(n)].add_vertex(v.id());
        }
    }
    std::cout << "3";
    auto new_size = (int)new_vertices.size();
    auto new_global_size = bulk::sum(new_world, new_size);
    std::cout << "4";
    auto new_H = pmondriaan::hypergraph(new_global_size, H.global_number_nets(),
                                        new_vertices, new_nets);
    std::cout << "5";
    remove_free_nets(new_world, new_H);
    std::cout << "6";

    return new_H;
}


} // namespace pmondriaan
