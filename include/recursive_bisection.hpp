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
                      long k,
                      double epsilon,
                      double eta,
                      pmondriaan::options opts,
                      std::string breaking_mode = "none",
                      std::string limit_edge_size = "false",
                      std::string simplify_mode = "complete");

std::vector<long>
compute_max_global_weight(long k_, long k_low, long k_high, long weight_mypart, long maxweight);

/**
 * Redistributes the hypergraph such that processors with my_part 0 contain all
 * vertices with label_low and all with my_part 1 contain all vertices with
 * label_high. Returns the end of part 0 if part 1 is not assigned any
 * processors.
 */
long redistribute_hypergraph(bulk::world& world,
                             pmondriaan::hypergraph& H,
                             long my_part,
                             long label_low,
                             long label_high,
                             long max_local_weight,
                             long weight_part_0,
                             long weight_part_1,
                             long p_low);

/**
 * Distributes the surplus vertices of a processor over the other processors.
 * Returns 0 if all surplus is gone, -1 if there is still surplus left.
 */
long reduce_surplus(bulk::world& world,
                    pmondriaan::hypergraph& H,
                    long label,
                    bulk::coarray<long>& surplus,
                    bulk::queue<long, long, long, long[]>& q,
                    bulk::queue<long, long>& cost_queue);

/**
 * Reorders the hypergraph such that all vertices with label_high are at the end of the vertex list.
 */
void reorder_hypergraph(pmondriaan::hypergraph& H, long start, long& end, long label_low, long label_high);

/**
 * Removes all cut nets from the hypergraph when the cut net metric is used.
 * The nets are stored in the cut_nets vector, to later add them to the hypergraph again.
 */
void remove_cut_nets(bulk::world& world,
                     bulk::world& sub_world,
                     pmondriaan::hypergraph& H,
                     std::vector<pmondriaan::net>& cut_nets);

/**
 * Adds all cut nets to the hypergraph again when the cut net metric is used.
 */
void add_cut_nets(bulk::world& world,
                  pmondriaan::hypergraph& H,
                  std::vector<pmondriaan::net>& cut_nets);

} // namespace pmondriaan
