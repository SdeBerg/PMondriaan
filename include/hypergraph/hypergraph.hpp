#pragma once

#include <iostream>
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
    vertex(int id, std::vector<int> nets, long weight = 1)
    : id_(id), nets_(nets), weight_(weight) {
        part_ = -1;
    }

    int id() { return id_; }
    std::vector<int>& nets() { return nets_; }
    long weight() { return weight_; }
    int part() { return part_; }
    auto degree() { return nets_.size(); }

    void set_id(int value) { id_ = value; }
    void add_weight(long value) { weight_ += value; }
    void set_part(int value) { part_ = value; }

    void remove_net(int n);

  private:
    int id_;
    std::vector<int> nets_;
    long weight_;
    int part_;
};

/**
 * A vertex has an id, a weight of type T and a list of the nets it is contained in.
 */
class net {
  public:
    net(int id, std::vector<int> vertices, long cost = 1)
    : id_(id), vertices_(vertices), cost_(cost) {
        global_size_ = 0;
    }

    int id() { return id_; }
    std::vector<int>& vertices() { return vertices_; }
    long cost() { return cost_; }
    auto size() { return vertices_.size(); }
    auto global_size() { return global_size_; }

    void set_global_size(size_t size) { global_size_ = size; }
    void add_vertex(int v) { vertices_.push_back(v); }

  private:
    int id_;
    std::vector<int> vertices_;
    long cost_;
    size_t global_size_;
};

/**
 * A hypergraph.
 */
class hypergraph {

  public:
    hypergraph(int global_size,
               std::vector<pmondriaan::vertex> vertices,
               std::vector<pmondriaan::net> nets)
    : global_size_(global_size), vertices_(std::move(vertices)),
      nets_(std::move(nets)) {
        global_to_local = std::unordered_map<int, int>();
        update_map();
    }

    // computes the total weight of the vertices
    long total_weight();

    // computes the sum of the weights of vertices in part
    long weight_part(int part);

    // computes the weights of all parts upto k
    std::vector<long> weight_all_parts(int k);

    // adds vertex to all nets
    void add_to_nets(pmondriaan::vertex& v);

    // removes vertex id from all nets
    void remove_from_nets(int id);

    // moves a vertex to the other part in 0,1
    void move(int id);

    // moves a vertex to the other part in 0,1 and adjusts the vector with counts of the parts
    void move(int id, std::vector<std::vector<long>>& C);

    // updates the global_to_local map
    void update_map();

    int local_id(int global_id) { return global_to_local[global_id]; }
    bool is_local(int global_id) {
        return (global_to_local.count(global_id) > 0);
    }

    void set_global_net_sizes(std::vector<size_t>& sizes);

    std::vector<pmondriaan::vertex>& vertices() { return vertices_; }
    pmondriaan::vertex& operator()(int index) { return vertices_[index]; }

    std::vector<pmondriaan::net>& nets() { return nets_; }
    pmondriaan::net& net(int index) { return nets_[index]; }

    auto size() { return vertices_.size(); }
    auto global_size() const { return global_size_; }
    auto& map() { return global_to_local; }

    void print();

  private:
    int global_size_;
    std::vector<pmondriaan::vertex> vertices_;
    std::vector<pmondriaan::net> nets_;
    std::unordered_map<int, int> global_to_local;
};

/**
 * Compute the global weight of a hypergraph.
 */
long global_weight(bulk::world& world, pmondriaan::hypergraph& H);

/**
 * Compute the global load imbalance of a hypergraph.
 */
double load_balance(bulk::world& world, pmondriaan::hypergraph& H, int k);

/**
 * Compute the cutsize with the correct metric of a local hypergraph
 */
long cutsize(pmondriaan::hypergraph& H, pmondriaan::m metric);

/**
 * Compute the cutsize with the correct metric of a distributed hypergraph
 */
long cutsize(bulk::world& world, pmondriaan::hypergraph& H, pmondriaan::m metric);

/**
 * Compute the global net sizes of a hypergraph.
 */
std::vector<size_t> global_net_sizes(bulk::world& world, pmondriaan::hypergraph& H);

/**
 * Removes all free nets.
 */
void remove_free_nets(bulk::world& world, pmondriaan::hypergraph& H);

/**
 * Creates a new hypergraph that only contains the vertices of H with local id between start and end.
 */
pmondriaan::hypergraph
create_new_hypergraph(bulk::world& world, pmondriaan::hypergraph& H, int start, int end);

/**
 * Creates a copy of a hypergraph and returns that copy.
 */
pmondriaan::hypergraph copy_hypergraph(pmondriaan::hypergraph& H);


} // namespace pmondriaan