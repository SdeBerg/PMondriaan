#pragma once

#include <vector>
#include <math.h>

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
		
		std::vector<int> nets() { return nets_; }
		
		long weight() { return weight_; }
		
		int part() { return part_; }
		
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
		
		void add_vertex(int v) { vertices_.push_back(v); }
		
		long cost() { return cost_; }
		
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
		hypergraph(std::vector<pmondriaan::vertex> vertices, std::vector<pmondriaan::net> nets)
			: vertices_(std::move(vertices)), nets_(std::move(nets)) {}
			
	long total_weight() {
		long total = 0;
		for (auto v : vertices_) {
			total += v.weight();
		}
		return total;
	}
			
		pmondriaan::vertex operator()(int index) {
			return vertices_[index];
		}
			
		auto size() const { return vertices_.size(); }

	private:
		std::vector<pmondriaan::vertex> vertices_;
		std::vector<pmondriaan::net> nets_;
};
	
	
} // namespace pmondriaan