#pragma once

#include <vector>

#include "hypergraph/hypergraph.hpp"
#include "options.hpp"

namespace pmondriaan {

/**
 * Returns a vector of ns randomly selected sample vertices.
 */
std::vector<int> sample_random(pmondriaan::hypergraph& H, int ns);

/**
 * Returns a vector of ns samples seleccted using the label propagation algorithm.
 */
std::vector<int> sample_lp(pmondriaan::hypergraph& H, pmondriaan::options& opts);

} // namespace pmondriaan
