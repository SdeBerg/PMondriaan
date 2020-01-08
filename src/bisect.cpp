#include <stdlib.h>
#include <vector>
#include <math.h>

#include "bisect.hpp"
#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

/**
 * Randomly bisects a hypergraph under the balance constraint and returns the weights of the two parts.
 */
std::vector<long> bisect_random(bulk::world& world, pmondriaan::hypergraph& H, double epsilon, int k, int label_0, int label_1) {
	
	auto max_weight_parts = std::vector<long>(2);
	max_weight_parts[0] = ceil(((double)H.total_weight() * ceil((double)k/2.0) / (double)k) * (1.0 + epsilon));
	max_weight_parts[1] = ceil(((double)H.total_weight() * (double)(k/2) / (double)k) * (1.0 + epsilon));
	
	world.log("max weight 0: %d", max_weight_parts[0]);
	world.log("max weight 1: %d", max_weight_parts[1]);

	auto weight_parts = std::vector<long>(2);
	
	for (auto i = 0u; i < H.size(); i++) {
		int random = rand();
		int part = random%2;
		// if the max weight is exceeded, we add the vertex to the other part
		if (weight_parts[part] + H(i).weight() > max_weight_parts[part]) {
			part = (part + 1)%2;
		}
		
		if (part == 0) {
			H(i).set_part(label_0);
		}
		else {
			H(i).set_part(label_1);
		}
		weight_parts[part] += H(i).weight();
	}
	
	return weight_parts;
}

} // namespace pmondriaan
