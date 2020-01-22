#pragma once

#include <vector>

#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

/**
 * Returns a vector of s randomly selected sample vertices.
 */
std::vector<int> sample_random(pmondriaan::hypergraph& H, int s);

} // namespace pmondriaan
