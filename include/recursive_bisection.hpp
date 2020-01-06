#pragma once

#include <stdlib.h>
#include <vector>
#include <string>

#include "hypergraph.hpp"
#include "bisect.hpp"

namespace pmondriaan {

/**
 * Recursively bisects a hypergraph into k parts.
 */
void recursive_bisect(bulk::world& world, pmondriaan::hypergraph& H, std::string mode, int k, double epsilon) {
	
	int s = world.rank();
	int p = world.active_processors();
	// the processors working with processor s on the same part are procs_mypart[0],..., procs_mypart[1] - 1.
	auto procs_mypart = std::vector<int>(2);
	procs_mypart[0] = 0;
	procs_mypart[1] = p;
	// we need to make sure every processor uses a different seed to generate random numbers
	srand(s + 1);
	int label_low = 0;
	int label_high = k - 1;
	
	auto weight_parts = std::vector<long>(2);
	// while we need to give more than one label, we bisect the hypergraph
	while (label_high - label_low > 1) {
		if (mode == "random") {
			weight_parts = bisect_random(world, H, epsilon, label_high - label_low + 1);
		}
		
		//TODO: reorder the hypergraph
		label_high = 
		label_low = 
	}
	
	
	
	
}

} // namespace pmondriaan
