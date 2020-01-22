#include <stdlib.h>
#include <vector>
#include <math.h>

#include "bisect.hpp"
#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/sample.hpp"

namespace pmondriaan {

/**
 * Randomly bisects a hypergraph under the balance constraint and returns the weights of the two parts.
 */
std::vector<long> bisect_random(bulk::world& world, pmondriaan::hypergraph& H, long max_weight_0, long max_weight_1, 
				int start, int end, int label_0, int label_1) {
	
	auto max_weight_parts = std::vector<long>{max_weight_0, max_weight_1};

	auto weight_parts = std::vector<long>(2);
	
	for (auto i = start; i < end; i++) {
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

/**
 * Bisects a hypergraph using the multilevel framework.
 
std::vector<long> bisect_multilevel(bulk::world&, pmondriaan::hypergraph& H, long max_weight_0, long max_weight_1,
			int start, int end, int label_0 = 0, int label_1 = 1) {
	
	long nc = 0;

	while ((HC.global_size() > coarsening_nrvertices) && (nc < coarsening_maxrounds)) {
		//we first select s samples	
				
				
				
				
	}
				
				
}*/

} // namespace pmondriaan
