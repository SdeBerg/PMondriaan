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
 * TODO: make sure the correct world is used (so not always all procs).
 */
std::vector<long> bisect_multilevel(bulk::world& world, pmondriaan::hypergraph& H, long max_weight_0, long max_weight_1,
			int start, int end, int label_0, int label_1) {
	
	int s = world.rank();
	int p = world.active_processors();
	
	auto weight_parts = std::vector<long>(2);
	
	//long nc = 0;

	//while ((HC.global_size() > options.coarsening_nrvertices) && (nc < options.coarsening_maxrounds)) {
		//we first select ns samples	
		//auto indices_samples = sample_random(HC, options.sample_size);
		auto indices_samples = sample_random(H, 10);
		auto local_samples = std::vector<std::vector<int>>(indices_samples.size());
		for (auto i = 0u; i < indices_samples.size(); i++) {
			//local_samples[i] = HC(indices_samples[i]).nets();
			local_samples[i] = H(indices_samples[i]).nets();
		}
		
		//broadcast samples
		//auto samples = bulk::coarray<int[]>(world, p * options.sample_size);
		auto samples = bulk::coarray<int[]>(world, p * 10);
		for (auto t = 0; t < p; t++) {
			samples(t)[{s * p, (s + 1) * p}] = local_samples;
		}
		world.sync();
		
		if (s == 0) {
			for (auto t = 0; t < p * 10; t++)
			world.log("t: %d, samples[t][0]: %d", t, samples[t][0]);
		}
				
	//}
				
	return weight_parts;		
}

} // namespace pmondriaan
