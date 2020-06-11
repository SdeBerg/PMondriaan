#pragma once

#include <string>

#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

bool create_random_hypergraph(std::string file, long n, long m, long nz);

} // namespace pmondriaan
