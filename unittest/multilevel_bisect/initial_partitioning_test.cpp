#include "pmondriaan.hpp"

#include "gtest/gtest.h"

namespace pmondriaan {
namespace {

TEST(InitialPartitioning, initial_partitioning) {
    auto H =
    pmondriaan::read_hypergraph("../test/data/matrices/dolphins/dolphins.mtx", "one")
    .value();
    std::mt19937 rng(1);
    pmondriaan::options opts;
    opts.KLFM_max_passes = 1;
    opts.metric = pmondriaan::m::cut_net;

    auto sol = pmondriaan::initial_partitioning(H, 31, 31, opts, rng);
    ASSERT_LE(H.weight_part(0), 31);
    ASSERT_LE(H.weight_part(1), 31);
    ASSERT_GE(sol, 6); 
}

} // namespace
} // namespace pmondriaan
