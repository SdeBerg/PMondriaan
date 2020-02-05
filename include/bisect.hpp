#pragma once

#include <vector>
#include <string>

#include "hypergraph/hypergraph.hpp"
#include "options.hpp"

namespace pmondriaan {

/**
 * Randomly bisects a hypergraph under the balance constraint and returns the weights of the two parts.
 */
std::vector<long> bisect_random(bulk::world& world, pmondriaan::hypergraph& H, long max_weight_0, long max_weight_1, 
			int start, int end, int label_0 = 0, int label_1 = 1);
			
/**
 * Bisects a hypergraph using the multilevel framework.
 */
std::vector<long> bisect_multilevel(bulk::world& world, pmondriaan::hypergraph& H, pmondriaan::options& opts, std::string sampling_mode,
			std::string metric, long max_weight_0, long max_weight_1, int start, int end, int label_0 = 0, int label_1 = 1);

} // namespace pmondriaan