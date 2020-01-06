#pragma once

#include <stdlib.h>
#include <vector>

#include "hypergraph.hpp"

namespace pmondriaan {

/**
 * Randomly bisects a hypergraph under the balance constraint and returns the weights of the two parts.
 */
std::vector<long> bisect_random(bulk::world& world, pmondriaan::hypergraph& H, double epsilon, int k = 2) {
	
	double max_weight = (double)H.total_weight() * math.ceil((double)k/2.0) / (double)k;
	long max_weight_part = avg_weight * (1.0 + epsilon);
	
	if (max_weight_part * 2.0 < avg_weight * 2.0) {
		max_weight_part++;
	}
	world.log("%d", max_weight_part);

	auto weight_parts = std::vector<long>(2);
	
	for (auto i = 0u; i < H.size(); i++) {
		int random = rand();
		int part = random%2;
		// if the max weight is exceeded, we add the vertex to the other part
		if (weight_parts[part] + H(i).weight() > max_weight_part) {
			part = (part + 1)%2;
		}
		
		H(i).set_part(part);
		weight_parts[part] += H(i).weight();
	}
	
	return weight_parts;
}

} // namespace pmondriaan
