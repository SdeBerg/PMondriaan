#include <stdlib.h>
#include <vector>
#include <string>

#include "recursive_bisection.hpp"
#include "hypergraph/hypergraph.hpp"
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
	
	// part 0 will always have the heaviest weight
	bulk::var<long> weight_part_0(world);
	bulk::var<long> weight_part_1(world);
	// while we need to give more than one label, we bisect the hypergraph
	while (label_high - label_low > 1) {
		int k_ = label_high - label_low + 1;
		if (mode == "random") {
			auto weight_parts = bisect_random(world, H, epsilon, k_, label_low, label_high);
			weight_part_0 = weight_parts[0];
			weight_part_1 = weight_parts[1];
		}
				
		world.sync();
		
		//TODO: implement this foldl version
		auto total_weight_part_0 = bulk::foldl(weight_part_0, procs_mypart[0], procs_mypart[1]);
		auto total_weight_part_1 = bulk::foldl(weight_part_1, procs_mypart[0], procs_mypart[1]);
		
		

		//TODO: reorder the hypergraph
		label_high = 
		label_low = 
	}
	
	
	
	
}

} // namespace pmondriaan
