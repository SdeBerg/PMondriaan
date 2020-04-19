#pragma once

#include <random>
#include <vector>

#include "hypergraph/hypergraph.hpp"
#include "options.hpp"

namespace pmondriaan {

/**
 * Returns a vector of ns randomly selected sample vertices.
 */
std::vector<long> sample_random(pmondriaan::hypergraph& H, long ns, std::mt19937& rng);

/**
 * Returns a vector of ns samples seleccted using the label propagation algorithm.
 */
std::vector<long>
sample_lp(pmondriaan::hypergraph& H, pmondriaan::options& opts, std::mt19937& rng);

} // namespace pmondriaan
