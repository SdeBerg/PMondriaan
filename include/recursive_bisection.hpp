#pragma once

#include <stdlib.h>
#include <vector>
#include <string>

#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

/**
 * Recursively bisects a hypergraph into k parts.
 */
void recursive_bisect(bulk::world& world, pmondriaan::hypergraph& H, std::string mode, int k, double epsilon, double eta);

/**
 * Redistributes the hypergraph such that processors with my_part 0 contain all vertices with label_low and
 * all with my_part 1 contain all vertices with label_high. Returns the end of part 0 if part 1 is not assigned
 * any processors.
 */
int redistribute_hypergraph(bulk::world& world, pmondriaan::hypergraph& H, std::vector<int> procs_mypart, int my_part, 
				int label_low, int label_high, long max_local_weight, long weight_part_0, long weight_part_1, int p_low);

/**
 * Distributes the surplus vertices of a processor over the other processors. Returns 0 if
 * all surplus is gone, -1 if there is still surplus left.
 */								
int reduce_surplus(bulk::world& world, pmondriaan::hypergraph& H, std::vector<int> procs_mypart, 
				int label, bulk::coarray<long>& surplus, bulk::queue<int, long, int, int[]>& q);

} // namespace pmondriaan
