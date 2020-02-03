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

long hypergraph::total_weight() {
	long total = 0;
	for (auto v : vertices_) {
		total += v.weight();
	}
	return total;
}

long hypergraph::global_weight(bulk::world& world) {
	
	bulk::var<long> local_weight(world);
	local_weight = total_weight();
	
	long global_weight = bulk::foldl(local_weight, [](auto& lhs, auto rhs) { lhs += rhs; });
	return global_weight;
}

//computes the sum of the weights of vertices in part
long hypergraph::weight_part(int part) {
	long total = 0;
	for (auto v : vertices_) {
		if(v.part() == part) {
			total += v.weight();
		}
	}
	return total;
}

//computes the weights of all parts upto kbhit
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

//removes id from all nets
void hypergraph::remove_from_nets(int id) {
	for (auto& net : nets_) {
		auto it = std::find(net.vertices().begin(), net.vertices().end(), id);
		if (it != net.vertices().end()) {
			std::iter_swap(it, net.vertices().end() - 1);
			net.vertices().erase(net.vertices().end() - 1);
		}
	}
}

void hypergraph::renumber_vertices() {
	for (auto i = 0u; i < vertices_.size(); i++) {
		global_to_local[vertices_[i].id()] = i;
	}
}

/**
 * Compute the global load imbalance of a hypergraph split into k parts.
 */
double compute_load_balance(bulk::world& world, pmondriaan::hypergraph& H, int k) {
	
	auto weight_parts_var = bulk::coarray<long>(world, k);
	auto weight_parts = H.weight_all_parts(k);
	
	for (int i = 0; i < k; i++) {
		weight_parts_var[i] = weight_parts[i];
	}
	
	world.sync();
	
	long global_weight = H.global_weight(world);
	
	//compute the global part weights
	weight_parts = pmondriaan::foldl(weight_parts_var, [](auto& lhs, auto rhs) { lhs += rhs; });
	
	//we compute the global part with largest weight
	long max_weight_part = *std::max_element(weight_parts.begin(), weight_parts.end());
	
	double eps = ((double)(max_weight_part*k) / (double)global_weight) - 1.0;
	
	return eps;
}

/**
 * Creates a new hypergraph that only contains the vertices of H with local id between start and end.
 */
pmondriaan::hypergraph create_new_hypergraph(pmondriaan::hypergraph& H, int start, int end) {
	
	std::vector<pmondriaan::vertex> new_vertices(H.vertices().begin() + start, H.vertices().begin() + end);
	auto new_nets = std::vector<pmondriaan::net>();
	for (auto& n : H.nets()) {
		new_nets.push_back(pmondriaan::net(n.id(), std::vector<int>()));
	}
	for (auto& v : new_vertices) {
		for (auto n : v.nets()) {
			new_nets[n].add_vertex(v.id());
		}
	}
	
	return pmondriaan::hypergraph(end - start, new_vertices, new_nets);
}
	
	
} // namespace pmondriaan