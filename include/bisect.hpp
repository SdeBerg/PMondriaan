#pragma once

#include <vector>

#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

/**
 * Randomly bisects a hypergraph under the balance constraint and returns the weights of the two parts.
 */
std::vector<long> bisect_random(bulk::world& world, pmondriaan::hypergraph& H, double epsilon, int k = 2, int label_0 = 0, int label_1 = 1);

} // namespace pmondriaan
