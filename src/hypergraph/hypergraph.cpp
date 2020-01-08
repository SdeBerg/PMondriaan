#include <vector>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

long hypergraph::total_weight() {
	long total = 0;
	for (auto v : vertices_) {
		total += v.weight();
	}
	return total;
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
	bulk::var<long> local_weight(world);
	local_weight = H.total_weight();
	
	auto weight_parts_var = std::vector<bulk::var<long>>(k);
	auto weight_parts = H.weight_all_parts(k);
	
	for (int i = 0; i < k; i++) {
		weight_parts_var[i] = weight_parts[i];
	}
	
	world.sync();
	
	long global_weight = bulk::foldl(local_weight, [](auto& lhs, auto rhs) { lhs += rhs; });
	
	//we compute the global part with largest weight
	long max_weight_part = bulk::foldl(weight_parts[0], [](auto& lhs, auto rhs) { lhs += rhs; });
	for (int i = 1; i < k; i++) {
		long w = bulk::foldl(weight_parts[i], [](auto& lhs, auto rhs) { lhs += rhs; } );
		if (w > max_weight_part) {
			max_weight_part = w;
		}
	}
	
	double eps = (double)max_weight_part / (double)global_weight * (double)k - 1.0;
	return eps;
}
	
	
} // namespace pmondriaan