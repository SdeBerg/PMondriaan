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
	
	
} // namespace pmondriaan