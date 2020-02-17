#include <stdlib.h>
#include <vector>
#include <math.h>

#include "bisect.hpp"
#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/coarsen.hpp"
#include "options.hpp"
#include "multilevel_bisect/initial_partitioning.hpp"

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
 * TODO: make sure the correct world is used (so not always all procs).
 */
std::vector<long> bisect_multilevel(bulk::world& world, pmondriaan::hypergraph& H, pmondriaan::options& opts, std::string sampling_mode,
			std::string metric, long max_weight_0, long max_weight_1, int start, int end, int label_0, int label_1) {
	
	auto weight_parts = std::vector<long>(2);
	
	auto H_reduced = pmondriaan::create_new_hypergraph(world, H, start, end);
	
	//int s = world.rank();
	//int p = world.active_processors();
	
	long nc = 0;
	
	auto HC_list = std::vector<pmondriaan::hypergraph>();
	HC_list.push_back(H_reduced);
	
	while ((HC_list[nc].global_size() > opts.coarsening_nrvertices) && (nc < opts.coarsening_maxrounds)) {
		HC_list.push_back(coarsen_hypergraph(world, HC_list[nc], opts, sampling_mode));
		nc++;
		world.log("After iteration %d, size is %d", nc, HC_list[nc].global_size());
	}
	
	pmondriaan::initial_partitioning(world, HC_list[nc], max_weight_0, max_weight_1, label_0, label_1);
				
	return weight_parts;		
}

} // namespace pmondriaan
