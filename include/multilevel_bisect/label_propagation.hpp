#include <vector>

#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

/**
 * Performs label propagation to create l groups of labels on the hypergraph H and returns 
 * a vector with all labels.
 */
std::vector<int> label_propagation(pmondriaan::hypergraph& H, int l, int max_iter);

} // namespace pmondriaan
