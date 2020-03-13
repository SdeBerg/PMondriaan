#pragma once

#include <string>
#include <vector>
#include <random>

#include "hypergraph/hypergraph.hpp"
#include "util/interval.hpp"
#include "options.hpp"

namespace pmondriaan {

/**
 * Bisect a hypergraph using the given bisection method and returns the weights of the two parts.
 */
std::vector<long> bisect(bulk::world& world,
                         pmondriaan::hypergraph& H,
                         std::string bisect_mode,
                         std::string sampling_mode,
                         pmondriaan::options& opts,
                         std::string metric,
                         long max_weight_0,
                         long max_weight_1,
                         int start,
                         int end,
                         interval labels,
                         std::mt19937& rng);

/**
 * Randomly bisects a hypergraph under the balance constraint and returns the weights of the two parts.
 */
std::vector<long> bisect_random(bulk::world& world,
                                pmondriaan::hypergraph& H,
                                long max_weight_0,
                                long max_weight_1,
                                int start,
                                int end,
                                interval labels,
                                std::mt19937& rng);

/**
 * Bisects a hypergraph using the multilevel framework.
 */
std::vector<long> bisect_multilevel(bulk::world& world,
                                    pmondriaan::hypergraph& H,
                                    pmondriaan::options& opts,
                                    std::string sampling_mode,
                                    std::string metric,
                                    long max_weight_0,
                                    long max_weight_1,
                                    int start,
                                    int end,
                                    interval labels,
                                    std::mt19937& rng);

} // namespace pmondriaan