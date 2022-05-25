#include <algorithm>
#include <iostream>
#include <random>
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

#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>

namespace pmondriaan {

// Naive and inefficient way to hash nets for mapping.
struct VectorHasher {
    std::size_t operator()(const std::vector<long> &V) const {
        return boost::hash_range(V.begin(), V.end());
    }
};

void vertex::remove_net(long n) {
    auto it = std::find(nets_.begin(), nets_.end(), n);
    if (it != nets_.end()) {
        std::iter_swap(it, nets_.end() - 1);
        nets_.pop_back();
    }
}

void vertex::add_net(long n) {
    nets_.push_back(n);
}

void net::remove_vertex(long v) {
    auto it = std::find(vertices_.begin(), vertices_.end(), v);
    if (it != vertices_.end()) {
        std::iter_swap(it, vertices_.end() - 1);
        vertices_.pop_back();
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
long hypergraph::weight_part(long part) {
    long total = 0;
    for (auto v : vertices_) {
        if (v.part() == part) {
            total += v.weight();
        }
    }
    return total;
}

// computes the weights of all parts upto k
std::vector<long> hypergraph::weight_all_parts(long k) {
    auto total = std::vector<long>(k);
    for (auto v : vertices_) {
        total[v.part()] += v.weight();
    }
    return total;
}

// add a vertex
void hypergraph::add_vertex(long id, std::vector<long> nets, long weight) {
    vertices_.push_back(pmondriaan::vertex(id, nets, weight));
    global_to_local[id] = (long)vertices_.size() - 1;
}

// add a net if it does not exist yet
void hypergraph::add_net(long id, std::vector<long> vertices, long cost) {
    if (net_global_to_local.count(id) == 0) {
        nets_.push_back(pmondriaan::net(id, vertices, cost));
        net_global_to_local[id] = (long)nets_.size() - 1;
        for (auto v : vertices) {
            vertices_[global_to_local[v]].add_net(id);
        }
    }
}

void hypergraph::add_to_nets(pmondriaan::vertex& v) {
    for (auto net_id : v.nets()) {
        nets_[local_id_net(net_id)].vertices().push_back(v.id());
    }
}

// removes id from all nets
void hypergraph::remove_from_nets(long id) {
    for (auto& net : nets_) {
        auto it = std::find(net.vertices().begin(), net.vertices().end(), id);
        if (it != net.vertices().end()) {
            std::iter_swap(it, net.vertices().end() - 1);
            net.vertices().pop_back();
        }
    }
}

// removes a free vertex from the vertex list
void hypergraph::remove_free_vertex(long id) {
    auto index = global_to_local[id];
    std::iter_swap(vertices_.begin() + index, vertices_.end() - 1);
    vertices().pop_back();
    global_to_local.insert_or_assign(vertices_[index].id(), index);
    global_to_local.erase(id);
}

// removes a net and the net from all net lists of vertices
void hypergraph::remove_net_by_index(long index) {
    long id = nets_[index].id();
    for (auto& v : nets_[index].vertices()) {
        vertices_[local_id(v)].remove_net(id);
    }
    std::iter_swap(nets_.begin() + index, nets_.end() - 1);
    nets_.pop_back();
    net_global_to_local.insert_or_assign(nets_[index].id(), index);
    net_global_to_local.erase(id);
}

// sorts the vertices in the nets of part 0 and 1 on their part
void hypergraph::sort_vertices_on_part(std::vector<std::vector<long>>& C) {
    for (auto i = 0u; i < nets_.size(); i++) {
        auto& net = nets_[i];
        auto& vertex_list = net.vertices();
        long index = 0;
        long end = net.size() - 1;
        while (index < C[i][0]) {
            if (vertices_[local_id(vertex_list[index])].part() != 0) {
                while (vertices_[local_id(vertex_list[end])].part() == 1) {
                    end--;
                }
                std::iter_swap(vertex_list.begin() + index, vertex_list.begin() + end);
                end--;
            } else {
                index++;
            }
        }
    }
}

// moves a vertex to the other part in 0,1 (before updating the counts)
void hypergraph::move_sorted(long id, std::vector<std::vector<long>>& C) {
    long idl = this->local_id(id);
    auto& vertex = vertices_[idl];
    vertex.set_part((vertex.part() + 1) % 2);
    if (vertex.part() == 0) {
        for (auto i = 0u; i < vertex.nets().size(); i++) {
            auto n = vertex.nets()[i];
            auto& vertex_list = this->net(n).vertices();
            long index = C[local_id_net(n)][0];
            while (vertex_list[index] != id) {
                index++;
            }
            std::iter_swap(vertex_list.begin() + index,
                           vertex_list.begin() + C[local_id_net(n)][0]);
        }
    } else {
        for (auto i = 0u; i < vertex.nets().size(); i++) {
            auto n = vertex.nets()[i];
            auto& vertex_list = this->net(n).vertices();
            long index = 0;
            while (vertex_list[index] != id) {
                index++;
            }
            std::iter_swap(vertex_list.begin() + index,
                           vertex_list.begin() + C[local_id_net(n)][0] - 1);
        }
    }
}

void hypergraph::move(long id, std::vector<std::vector<long>>& C) {
    auto& v = vertices_[this->local_id(id)];
    long from = v.part();
    move_sorted(id, C);
    long to = v.part();
    for (auto n : v.nets()) {
        C[local_id_net(n)][from]--;
        C[local_id_net(n)][to]++;
    }
}

void hypergraph::move(long id,
                      std::vector<std::vector<long>>& C,
                      std::vector<std::vector<long>>& C_loc) {
    auto& v = vertices_[this->local_id(id)];
    long from = v.part();
    move_sorted(id, C_loc);
    long to = v.part();
    for (auto n : v.nets()) {
        C[local_id_net(n)][from]--;
        C[local_id_net(n)][to]++;
        C_loc[local_id_net(n)][from]--;
        C_loc[local_id_net(n)][to]++;
    }
}

void hypergraph::update_map() {
    global_to_local.clear();
    for (auto i = 0u; i < vertices_.size(); i++) {
        global_to_local[vertices_[i].id()] = i;
    }
}

void hypergraph::update_map_nets() {
    net_global_to_local.clear();
    for (auto i = 0u; i < nets_.size(); i++) {
        net_global_to_local[nets_[i].id()] = i;
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

// For testing purposes, checks if the maps are correct
void hypergraph::check_maps() {
    for (auto v : vertices_) {
        if (global_to_local.find(v.id()) == global_to_local.end()) {
            std::cout << "Vertex not found!!\n";
        } else if (vertices_[local_id(v.id())].id() != v.id()) {
            std::cout << "Local id vertex incorrect!\n";
        }
        for (auto n : v.nets()) {
            if (net_global_to_local.find(n) == net_global_to_local.end()) {
                std::cout << "Net not found!\n";
            } else if (nets_[local_id_net(n)].id() != n) {
                std::cout << "Local id net incorrect!\n";
            }
        }
    }
    for (auto n : nets_) {
        for (auto v : n.vertices()) {
            if (global_to_local.find(v) == global_to_local.end()) {
                std::cout << "Vertex not found!!\n";
            } else if (vertices_[local_id(v)].id() != v) {
                std::cout << "Local id vertex incorrect!\n";
            }
        }
        if (net_global_to_local.find(n.id()) == net_global_to_local.end()) {
            std::cout << "Net not found!\n";
        } else if (nets_[local_id_net(n.id())].id() != n.id()) {
            std::cout << "Local id net incorrect!\n";
        }
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
 * Initialize the counts for parts 0,1 for a parallel hypergraph.
 */
std::vector<std::vector<long>> init_counts(bulk::world& world, pmondriaan::hypergraph& H) {
    auto s = world.rank();

    auto local_counts = init_counts(H);
    auto net_partition =
    bulk::block_partitioning<1>({H.global_number_nets()},
                                {(size_t)world.active_processors()});
    auto count_queue = bulk::queue<long, long>(world);
    // We send all counts of part 0 that are greater than 0 to the responsible processor
    for (auto i = 0u; i < local_counts.size(); i++) {
        if (local_counts[i][0] > 0) {
            auto global_id = H.global_id_net(i);
            count_queue(net_partition.owner(global_id)).send(global_id, local_counts[i][0]);
        }
    }
    world.sync();

    auto C_my_nets = bulk::coarray<long>(world, net_partition.local_count(s));
    for (auto i = 0u; i < net_partition.local_count(s); i++) {
        C_my_nets[i] = 0;
    }
    for (const auto& [net, count] : count_queue) {
        C_my_nets[net_partition.local(net)[0]] += count;
    }

    // The following comment contains the more elegant version that seems to be inefficient on some versions of MPI.

    /*// We get the correct values for the needed counts
    auto future_counts = std::vector<bulk::future<long>>();
    for (auto i = 0u; i < local_counts.size(); i++) {
        auto net = H.nets()[i].id();
        future_counts.push_back(
        C_my_nets(net_partition.owner(net))[net_partition.local(net)[0]].get());
    }

    world.sync();

    for (auto i = 0u; i < local_counts.size(); i++) {
        local_counts[i][0] = future_counts[i].value();
        local_counts[i][1] = (long)H.nets()[i].global_size() - future_counts[i].value();
    }*/

    auto request_queue = bulk::queue<int, int, long>(world);
    for (auto i = 0u; i < local_counts.size(); i++) {
        auto net = H.nets()[i].id();
        request_queue(net_partition.owner(net)).send(s, i, net);
    }
    world.sync();

    auto response_queue = bulk::queue<int, long>(world);

    for (const auto& [sender, idx, net] : request_queue) {
        response_queue(sender).send(idx, C_my_nets[net_partition.local(net)[0]]);
    }

    world.sync();

    for (const auto& [i, countresponse] : response_queue) {
        local_counts[i][0] = countresponse;
        local_counts[i][1] = (long)H.nets()[i].global_size() - countresponse;
    }

    return local_counts;
}

/**
 * Recompute the global size of a hypergraph.
 */
void recompute_global_size(bulk::world& world, pmondriaan::hypergraph& H) {
    auto new_size = bulk::sum(world, H.size());
    H.set_global_size(new_size);
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
 * Compute the global weight of a part of a hypergraph.
 */
long global_weight_part(bulk::world& world, pmondriaan::hypergraph& H, int part) {
    auto local_weight = H.weight_part(part);
    long global_weight = bulk::sum(world, local_weight);
    return global_weight;
}

std::vector<long> global_weight_parts(bulk::world& world, pmondriaan::hypergraph& H, long k) {
    auto weight_parts_coar = bulk::coarray<long>(world, k);
    auto weight_parts = H.weight_all_parts(k);
    for (long i = 0; i < k; i++) {
        weight_parts_coar[i] = weight_parts[i];
    }
    weight_parts =
    bulk::foldl_each(weight_parts_coar, [](auto& lhs, auto rhs) { lhs += rhs; });
    return weight_parts;
}

/**
 * Compute the load imbalance of a hypergraph.
 */
double load_balance(pmondriaan::hypergraph& H, long k) {
    auto weight_parts = H.weight_all_parts(k);
    // we compute the global part with largest weight
    long max_weight_part = *std::max_element(weight_parts.begin(), weight_parts.end());
    double eps = ((double)(max_weight_part * k) / (double)H.total_weight()) - 1.0;
    return eps;
}

/**
 * Compute the global load imbalance of a hypergraph split into k parts.
 */
double load_balance(bulk::world& world, pmondriaan::hypergraph& H, long k) {
    auto weight_parts = global_weight_parts(world, H, k);
    auto global_weight = pmondriaan::global_weight(world, H);

    long sum_weight_parts = 0;
    for (auto w : weight_parts) {
        sum_weight_parts += w;
    }
    if (sum_weight_parts != global_weight) {
        std::cerr << "Global weight not equal to sum of part weights!";
    }

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
            auto labels_net = std::unordered_set<long>();
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
            auto labels_net = std::unordered_set<long>();
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
    auto net_partition =
    bulk::block_partitioning<1>({H.global_number_nets()},
                                {(size_t)world.active_processors()});

    // this queue contains all labels present for each net and its cost
    auto labels = bulk::queue<long, long[], long>(world);
    for (auto& net : H.nets()) {
        auto labels_net = std::unordered_set<long>();
        for (auto& v : net.vertices()) {
            labels_net.insert(H(H.local_id(v)).part());
        }
        labels(net_partition.owner(net.id()))
        .send(net.id(), std::vector<long>(labels_net.begin(), labels_net.end()),
              net.cost());
    }
    world.sync();

    auto total_cut =
    std::vector<std::unordered_set<long>>(net_partition.local_size(world.rank())[0]);
    auto cost_nets = std::vector<long>(net_partition.local_size(world.rank())[0], 0);
    for (const auto& [net, labels_net, cost] : labels) {
        cost_nets[net_partition.local({(size_t)net})[0]] = cost;
        for (auto l : labels_net) {
            total_cut[net_partition.local({(size_t)net})[0]].insert(l);
        }
    }

    switch (metric) {
    case pmondriaan::m::cut_net: {
        for (auto i = 0u; i < total_cut.size(); i++) {
            if (total_cut[i].size() > 1) {
                result += cost_nets[i];
            }
        }
        break;
    }
    case pmondriaan::m::lambda_minus_one: {
        for (auto i = 0u; i < total_cut.size(); i++) {
            if (total_cut[i].size() > 1) {
                result += (total_cut[i].size() - 1) * cost_nets[i];
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
 * Compute the cutsize of a bisected hypergraph using the vector C of the counts of all nets
 */
long cutsize(pmondriaan::hypergraph& H, std::vector<std::vector<long>>& C) {
    long cut = 0;
    for (auto i = 0u; i < C.size(); i++) {
        if ((C[i][0] > 0) && (C[i][1] > 0)) {
            cut += H.nets()[i].cost();
        }
    }
    return cut;
}

/**
 * Compute the global net sizes of a hypergraph.
 */
std::vector<size_t> global_net_sizes(bulk::world& world, pmondriaan::hypergraph& H) {

    auto net_partition =
    bulk::block_partitioning<1>({H.global_number_nets()},
                                {(size_t)world.active_processors()});

    auto nets = H.nets();
    bulk::queue<long, size_t> net_size_queue(world);

    for (auto i = 0u; i < nets.size(); i++) {
        net_size_queue(net_partition.owner(nets[i].id()))
        .send(nets[i].id(), nets[i].size());
    }
    world.sync();

    auto size_nets = std::vector<size_t>(net_partition.local_count(world.rank()), 0);

    for (const auto& [id, size] : net_size_queue) {
        auto local_id = net_partition.local({(size_t)id})[0];
        size_nets[local_id] += size;
    }

    for (auto i = 0u; i < size_nets.size(); i++) {
        if (size_nets[i] > 0) {
            for (int t = 0; t < world.active_processors(); t++) {
                net_size_queue(t).send(net_partition.global(i, world.rank())[0],
                                       size_nets[i]);
            }
        }
    }
    world.sync();
    auto result = std::vector<size_t>(nets.size(), 0);
    for (const auto& [id, size] : net_size_queue) {
        if (H.is_local_net(id)) {
            result[H.local_id_net(id)] = size;
        }
    }
    H.set_global_net_sizes(result);
    return result;
}


/**
 * Removes all free nets.
 */
void remove_free_nets(bulk::world& world, pmondriaan::hypergraph& H, size_t max_size) {
    auto net_sizes = global_net_sizes(world, H);
    std::unordered_set<long> remove_nets;
    for (auto n = 0u; n < H.nets().size(); n++) {
        if (net_sizes[n] <= max_size || H.nets()[n].size() == 0) {
            remove_nets.insert(H.nets()[n].id());
        }
    }
    for (auto n : remove_nets) {
        H.remove_net_by_index(H.local_id_net(n));
    }
}

/**
 * Removes all free nets.
 */
void remove_free_nets(pmondriaan::hypergraph& H, size_t max_size) {
    std::unordered_set<long> remove_nets;
    for (auto n = 0u; n < H.nets().size(); n++) {
        if (H.nets()[n].size() <= max_size) {
            remove_nets.insert(H.nets()[n].id());
        } else {
            H.nets()[n].set_global_size(H.nets()[n].size());
        }
    }
    for (auto n : remove_nets) {
        H.remove_net_by_index(H.local_id_net(n));
    }
}

/**
 * Simplifies all duplicate nets.
 */
void simplify_duplicate_nets(pmondriaan::hypergraph& H) {
    // Initialize maps
    boost::unordered_map<std::vector<long>, long, VectorHasher> edge_counts;
    std::vector<long> remove_nets;

    // Simplify nets that are not too large.
    // Note that smaller hyperedges are more common in most cases, if not then simplification
    // is of questionable use anyways.
    // Other than this we just sum the costs for some set of vertices found to map to
    // a set of identifyer long integers in the pins of a hyperedges.
    for (auto& n : H.nets()) {
        if(n.vertices().size() < 100) {
            if(edge_counts.find(n.vertices()) == edge_counts.end()) {
                edge_counts[n.vertices()] = n.id();
            }
            else {
                auto& net = H.nets()[H.local_id_net(edge_counts[n.vertices()])];
                net.set_cost(n.cost() + net.cost());
                remove_nets.push_back(n.id());
            }
        }
    }

    for (auto n : remove_nets) {
        H.remove_net_by_index(H.local_id_net(n));
    }

}

/**
 * Break triples, this does not account for the policy around the usage..
 */
void break_triples(pmondriaan::hypergraph& H) {
    std::vector<long> remove_nets;

    long cid = 1;

    for (auto n = 0u; n < H.nets().size(); n++) {
        if (H.nets()[n].id() >= cid) {
            cid = H.nets()[n].id() + 1;
        }
    }

    // Needs to be done through indices because we add nets.
    for (auto n = 0u; n < H.nets().size(); n++) {
        long cost = H.nets()[n].cost();
        if(H.nets()[n].size() != 3) {
            H.nets()[n].set_cost(cost * 2);
        }
    }

    for (auto n = 0u; n < H.nets().size(); n++) {
        long cost = H.nets()[n].cost();
        if(H.nets()[n].size() == 3) {
            remove_nets.push_back(H.nets()[n].id());
            std::vector<long> verts1;
            std::vector<long> verts2;
            std::vector<long> verts3;
            verts1.push_back(H.nets()[n].vertices()[0]);
            verts1.push_back(H.nets()[n].vertices()[1]);
            verts2.push_back(H.nets()[n].vertices()[1]);
            verts2.push_back(H.nets()[n].vertices()[2]);
            verts3.push_back(H.nets()[n].vertices()[0]);
            verts3.push_back(H.nets()[n].vertices()[2]);
            H.add_net(cid + 0, verts1, cost);
            H.add_net(cid + 1, verts2, cost);
            H.add_net(cid + 2, verts3, cost);
            cid += 3;
        }
    }

    for (auto n : remove_nets) {
        H.remove_net_by_index(H.local_id_net(n));
    }

    simplify_duplicate_nets(H);
}

/**
 * Sorts vertices.
 */
void sort_verts(pmondriaan::hypergraph& H) {
    std::sort(H.vertices().begin(), H.vertices().end());
}

/**
 * Creates a new hypergraph that only contains the vertices of H with local id between start and end.
 */
pmondriaan::hypergraph
create_new_hypergraph(bulk::world& new_world, pmondriaan::hypergraph& H, long start, long end) {

    std::vector<pmondriaan::vertex> new_vertices(H.vertices().begin() + start,
                                                 H.vertices().begin() + end);
    auto new_nets = std::vector<pmondriaan::net>();
    for (auto& n : H.nets()) {
        new_nets.push_back(pmondriaan::net(n.id(), std::vector<long>()));
    }

    for (auto& v : new_vertices) {
        for (auto n : v.nets()) {
            new_nets[H.local_id_net(n)].add_vertex(v.id());
        }
    }

    auto new_size = new_vertices.size();
    auto new_global_size = bulk::sum(new_world, new_size);
    auto new_H = pmondriaan::hypergraph(new_global_size, H.global_number_nets(),
                                        new_vertices, new_nets);

    remove_free_nets(new_world, new_H, 1);
    return new_H;
}


} // namespace pmondriaan
