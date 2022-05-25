#pragma once

#include <random>
#include <string>
#include <vector>

#include "hypergraph/hypergraph.hpp"
#include "options.hpp"
#include "util/interval.hpp"

namespace pmondriaan {

/**
 * Bisect a hypergraph using the given bisection method and returns the weights of the two parts.
 */
std::vector<long> bisect(bulk::world& world,
                         pmondriaan::hypergraph& H,
                         pmondriaan::options& opts,
                         long max_weight_0,
                         long max_weight_1,
                         long start,
                         long end,
                         interval labels,
                         std::mt19937& rng,
                         std::string breaking_mode,
                         std::string limit_edge_size = "false",
                         std::string simplify_mode = "complete");

/**
 * Randomly bisects a hypergraph under the balance constraint and returns the weights of the two parts.
 */
std::vector<long> bisect_random(pmondriaan::hypergraph& H,
                                long max_weight_0,
                                long max_weight_1,
                                long start,
                                long end,
                                interval labels,
                                std::mt19937& rng);

/**
 * Bisects a hypergraph using the multilevel framework.
 */
std::vector<long> bisect_multilevel(bulk::world& world,
                                    pmondriaan::hypergraph& H,
                                    pmondriaan::options& opts,
                                    long max_weight_0,
                                    long max_weight_1,
                                    long start,
                                    long end,
                                    interval labels,
                                    std::mt19937& rng,
                                    std::string breaking_mode,
                                    std::string limit_edge_size = "false",
                                    std::string simplify_mode = "complete");

} // namespace pmondriaan
