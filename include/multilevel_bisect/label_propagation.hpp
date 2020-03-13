#include <random>
#include <vector>

#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

/**
 * Performs label propagation to create l groups of labels on the hypergraph H
 * and returns a vector with all labels.
 */
std::vector<int>
label_propagation(pmondriaan::hypergraph& H, int l, int max_iter, int min_size, std::mt19937& rng);

std::vector<int> label_propagation_bisect(pmondriaan::hypergraph& H,
                                          int max_iter,
                                          long max_weight_0,
                                          long max_weight_1,
                                          std::mt19937& rng);

} // namespace pmondriaan
