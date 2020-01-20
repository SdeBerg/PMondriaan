#pragma once

#include <vector>
#include <iostream>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

namespace pmondriaan {
	
/**
 * A vertex has an id, a weight of type T and a list of the nets it is contained in.
 */
class vertex {
	public:
		vertex(int id, std::vector<int> nets, long weight = 1) : id_(id), nets_(nets), weight_(weight) {}
		
		int id() { return id_; }
		std::vector<int>& nets() { return nets_; }
		long weight() { return weight_; }
		int part() { return part_; }
		
		void set_id(int value) {id_ = value; }
		void set_part(int value) { part_ = value; }
		
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
		net(int id, std::vector<int> vertices, long cost = 1) : id_(id), vertices_(vertices), cost_(cost) {}
		
		int id() { return id_; }
		std::vector<int> vertices() { return vertices_; }
		long cost() { return cost_; }
		
		void add_vertex(int v) { vertices_.push_back(v); }
		
	private:
		int id_;
		std::vector<int> vertices_;
		long cost_;
};

/**
 * A hypergraph.
 */
class hypergraph {
	
	public:
		hypergraph(int global_size, std::vector<pmondriaan::vertex> vertices, std::vector<pmondriaan::net> nets)
			: global_size_(global_size), vertices_(std::move(vertices)), nets_(std::move(nets)) {}
		
		//computes the total weight of the vertices
		long total_weight();
		
		//computes the sum of the weights of vertices in part
		long weight_part(int part);
		
		//computes the weights of all parts upto kbhit
		std::vector<long> weight_all_parts(int k);
			
		std::vector<pmondriaan::vertex>& vertices() { return vertices_; }
		pmondriaan::vertex& operator()(int index) { return vertices_[index]; }
		
		auto size() const { return vertices_.size(); }
		auto global_size() const {return global_size_; }

	private:
		int global_size_;
		std::vector<pmondriaan::vertex> vertices_;
		std::vector<pmondriaan::net> nets_;
};

/**
 * Compute the global load imbalance of a hypergraph.
 */
double compute_load_balance(bulk::world& world, pmondriaan::hypergraph& H, int k);
	
	
} // namespace pmondriaan