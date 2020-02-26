#pragma once

#include <stdlib.h>
#include <string>
#include <vector>

#include "hypergraph/hypergraph.hpp"
#include "options.hpp"
#include "util/interval.hpp"

namespace pmondriaan {

/**
 * Recursively bisects a hypergraph into k parts.
 */
void recursive_bisect(bulk::world& world,
                      pmondriaan::hypergraph& H,
                      std::string bisect_mode,
                      std::string sampling_mode,
                      std::string metric,
                      int k,
                      double epsilon,
                      double eta,
                      pmondriaan::options& opts);

std::vector<long>
compute_max_global_weight(int k_, int k_low, int k_high, long weight_mypart, long maxweight);

/**
 * Redistributes the hypergraph such that processors with my_part 0 contain all
 * vertices with label_low and all with my_part 1 contain all vertices with
 * label_high. Returns the end of part 0 if part 1 is not assigned any
 * processors.
 */
int redistribute_hypergraph(bulk::world& world,
                            pmondriaan::hypergraph& H,
                            std::vector<int> procs_mypart,
                            int my_part,
                            int label_low,
                            int label_high,
                            long max_local_weight,
                            long weight_part_0,
                            long weight_part_1,
                            int p_low);

/**
 * Distributes the surplus vertices of a processor over the other processors.
 * Returns 0 if all surplus is gone, -1 if there is still surplus left.
 */
int reduce_surplus(bulk::world& world,
                   pmondriaan::hypergraph& H,
                   std::vector<int> procs_mypart,
                   int label,
                   bulk::coarray<long>& surplus,
                   bulk::queue<int, long, int, int[]>& q);

/**
 * Reorders the hypergraph such that all vertices with label_high are at the end of the vertex list.
 */
void reorder_hypergraph(pmondriaan::hypergraph& H, int start, int& end, int label_low, int label_high);

} // namespace pmondriaan
