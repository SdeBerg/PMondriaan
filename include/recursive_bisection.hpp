#pragma once

#include <stdlib.h>
#include <vector>
#include <string>

#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

/**
 * Recursively bisects a hypergraph into k parts.
 */
void recursive_bisect(bulk::world& world, pmondriaan::hypergraph& H, std::string mode, int k, double epsilon);

} // namespace pmondriaan
