#pragma once

#include <cassert>
#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "options.hpp"
#include "util/interval.hpp"

namespace pmondriaan {

/**
 * A vertex has an id, a weight of type T and a list of the nets it is contained in.
 */
class vertex {
  public:
    vertex(long id, std::vector<long> nets, long weight = 1)
    : id_(id), nets_(nets), weight_(weight) {}

    long id() { return id_; }
    std::vector<long>& nets() { return nets_; }
    long weight() { return weight_; }
    long part() { return part_; }
    auto degree() { return nets_.size(); }

    void set_id(long value) { id_ = value; }
    void add_weight(long value) { weight_ += value; }
    void set_part(long value) { part_ = value; }

    void remove_net(long n);

  private:
    long id_;
    std::vector<long> nets_;
    long weight_;
    long part_ = -1;
};

/**
 * A vertex has an id, a weight of type T and a list of the nets it is contained in.
 */
class net {
  public:
    net(long id, std::vector<long> vertices, long cost = 1)
    : id_(id), vertices_(vertices), cost_(cost) {}

    long id() const { return id_; }

    std::vector<long>& vertices() { return vertices_; }
    const std::vector<long>& vertices() const { return vertices_; }

    long cost() const { return cost_; }
    auto size() const { return vertices_.size(); }
    auto global_size() const { return global_size_; }

    void set_global_size(size_t size) { global_size_ = size; }
    void add_vertex(long v) { vertices_.push_back(v); }

    double scaled_cost() const {
        return (double)cost_ / ((double)global_size_ - 1.0);
    }

  private:
    long id_;
    std::vector<long> vertices_;
    long cost_;
    size_t global_size_ = 0;
};

/**
 * A hypergraph.
 */
class hypergraph {

  public:
    hypergraph(size_t global_size,
               size_t global_number_nets,
               std::vector<pmondriaan::vertex> vertices,
               std::vector<pmondriaan::net> nets,
               size_t nr_nz = 0)
    : global_size_(global_size), global_number_nets_(global_number_nets),
      vertices_(std::move(vertices)), nets_(std::move(nets)), nr_nz_(nr_nz) {
        update_map();
        update_map_nets();
    }

    hypergraph(const hypergraph& other)
    : global_size_(other.global_size_), global_number_nets_(other.global_number_nets_),
      vertices_(other.vertices_), nr_nz_(other.nr_nz_) {
        for (const auto& n : other.nets()) {
            nets_.push_back(pmondriaan::net(n.id(), std::vector<long>()));
            nets_.back().set_global_size(n.global_size());
        }

        for (auto& v : vertices_) {
            for (auto n : v.nets()) {
                nets_[other.local_id_net(n)].add_vertex(v.id());
            }
        }

        update_map();
        update_map_nets();
    }

    hypergraph(hypergraph&& other) = default;

    // computes the total weight of the vertices
    long total_weight();

    // computes the sum of the weights of vertices in part
    long weight_part(long part);

    // computes the weights of all parts upto k
    std::vector<long> weight_all_parts(long k);

    // add a vertex
    void add_vertex(long id, std::vector<long> nets, long weight = 1);

    // add a net if it does not exist yet
    void add_net(long id, std::vector<long> vertices, long cost = 1);

    // adds vertex to all nets
    void add_to_nets(pmondriaan::vertex& v);

    // removes vertex id from all nets
    void remove_from_nets(long id);

    // removes a free vertex from the vertex list
    void remove_free_vertex(long id);

    // removes a net and the net from all net lists of vertices
    void remove_net_by_index(long index);

    // sorts the vertices in the nets of part 0 and 1 on their part
    void sort_vertices_on_part(std::vector<std::vector<long>>& C);

    // moves a vertex to the other part in 0,1
    void move_sorted(long id, std::vector<std::vector<long>>& C);

    // moves a vertex to the other part in 0,1 and adjusts the vector with counts of the parts
    void move(long id, std::vector<std::vector<long>>& C);

    // moves a vertex to the other part in 0,1 and adjusts the vector with counts of the parts for parallel hypergraph
    void move(long id, std::vector<std::vector<long>>& C, std::vector<std::vector<long>>& C_loc);

    // updates the global_to_local map
    void update_map();

    // updates the map for the nets
    void update_map_nets();

    long local_id(long global_id) {
        assert(global_to_local.find(global_id) != global_to_local.end());
        assert(global_to_local[global_id] >= 0 &&
               (size_t)global_to_local[global_id] < vertices_.size());
        return global_to_local[global_id];
    }

    long local_id_net(long global_id) const {
        return net_global_to_local.at(global_id);
    }
    long global_id_net(long local_id) { return nets_[local_id].id(); }
    bool is_local(long global_id) {
        return (global_to_local.count(global_id) > 0);
    }
    bool is_local_net(long global_id) {
        return net_global_to_local.count(global_id) > 0;
    }

    void set_global_size(size_t size) { global_size_ = size; }
    void set_global_net_sizes(std::vector<size_t>& sizes);

    std::vector<pmondriaan::vertex>& vertices() { return vertices_; }
    const std::vector<pmondriaan::vertex>& vertices() const {
        return vertices_;
    }
    pmondriaan::vertex& operator()(long index) {
        assert(index >= 0 && (size_t)index < vertices_.size());
        return vertices_[index];
    }

    std::vector<pmondriaan::net>& nets() { return nets_; }
    const std::vector<pmondriaan::net>& nets() const { return nets_; }

    pmondriaan::net& net(long id) {
        assert(net_global_to_local.find(id) != net_global_to_local.end());
        auto local_id = net_global_to_local[id];
        assert(local_id >= 0 && (size_t)local_id < nets_.size());
        return nets_[local_id];
    }

    auto size() { return vertices_.size(); }
    auto global_size() const { return global_size_; }
    auto global_number_nets() const { return global_number_nets_; }
    auto nr_nz() const { return nr_nz_; }
    auto& map() { return global_to_local; }
    auto& map_nets() { return net_global_to_local; }

    void print();
    // For testing purposes, checks if the maps are correct
    void check_maps();

  private:
    size_t global_size_;
    size_t global_number_nets_;
    std::vector<pmondriaan::vertex> vertices_;
    std::vector<pmondriaan::net> nets_;
    size_t nr_nz_;
    std::unordered_map<long, long> global_to_local;
    std::unordered_map<long, long> net_global_to_local;
};

/**
 * Initialize the counts for parts 0,1.
 */
std::vector<std::vector<long>> init_counts(pmondriaan::hypergraph& H);

/**
 * Initialize the counts for parts 0,1 for a parallel hypergraph.
 */
std::vector<std::vector<long>> init_counts(bulk::world& world, pmondriaan::hypergraph& H);

/**
 * Recompute the global size of a hypergraph.
 */
void recompute_global_size(bulk::world& world, pmondriaan::hypergraph& H);

/**
 * Compute the global weight of a hypergraph.
 */
long global_weight(bulk::world& world, pmondriaan::hypergraph& H);

/**
 * Compute the global weight of a part of a hypergraph.
 */
long global_weight_part(bulk::world& world, pmondriaan::hypergraph& H, int part);

/**
 * Compute the global weight of all parts up to k of a hypergraph.
 */
std::vector<long> global_weight_parts(bulk::world& world, pmondriaan::hypergraph& H, long k);

/**
 * Compute the load imbalance of a hypergraph.
 */
double load_balance(pmondriaan::hypergraph& H, long k);

/**
 * Compute the global load imbalance of a distributed hypergraph.
 */
double load_balance(bulk::world& world, pmondriaan::hypergraph& H, long k);

/**
 * Compute the cutsize with the correct metric of a local hypergraph
 */
long cutsize(pmondriaan::hypergraph& H, pmondriaan::m metric);

/**
 * Compute the cutsize with the correct metric of a distributed hypergraph
 */
long cutsize(bulk::world& world, pmondriaan::hypergraph& H, pmondriaan::m metric);

/**
 * Compute the cutsize of a bisected hypergraph using the vector C of the counts of all nets
 */
long cutsize(pmondriaan::hypergraph& H, std::vector<std::vector<long>>& C);

/**
 * Compute the global net sizes of a hypergraph.
 */
std::vector<size_t> global_net_sizes(bulk::world& world, pmondriaan::hypergraph& H);

/**
 * Removes all free nets.
 */
void remove_free_nets(bulk::world& world, pmondriaan::hypergraph& H, size_t max_size);

/**
 * Removes all free nets.
 */
void remove_free_nets(pmondriaan::hypergraph& H, size_t max_size);

/**
 * Creates a new hypergraph that only contains the vertices of H with local id between start and end.
 */
pmondriaan::hypergraph
create_new_hypergraph(bulk::world& world, pmondriaan::hypergraph& H, long start, long end);

} // namespace pmondriaan
