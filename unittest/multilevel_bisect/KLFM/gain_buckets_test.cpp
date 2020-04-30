#include "pmondriaan.hpp"

#include "gtest/gtest.h"

namespace pmondriaan {
namespace {

std::string mtx_three_nonzeros = R"(%%MatrixMarket matrix coordinate real general
3 3 4
1 1 1.0
2 1 1.0
2 2 1.0
1 3 1.0
)";

TEST(GainBucket, GainBucketInit) {
    std::stringstream mtx_ss(mtx_three_nonzeros);
    auto H = read_hypergraph_istream(mtx_ss, "one").value();
    H(0).set_part(0);
    H(1).set_part(1);
    H(2).set_part(0);
    auto C = std::vector<std::vector<long>>(2);
    C[0] = {2, 0};
    C[1] = {1, 1};
    std::mt19937 rng(1);
    auto g = pmondriaan::gain_structure(H, C);
    ASSERT_EQ(g.part_next(3, 3, rng), 1);
    ASSERT_EQ(g.gain_next(0), 0);
    ASSERT_EQ(g.gain_next(1), 1);
    ASSERT_EQ(g.next(1), 1);
    g.move(g.next(g.part_next(3, 3, rng)));
    ASSERT_EQ(g.part_next(3, 3, rng), 0);
    ASSERT_EQ(g.next(1), -1);
    ASSERT_EQ(g.next(0), 2);
    ASSERT_EQ(g.gain_next(0), -1);

    pmondriaan::options opts;
    opts.KLFM_max_passes = 1;
    opts.metric = pmondriaan::m::cut_net;
    auto sol = pmondriaan::initial_partitioning(H, 3, 3, opts, rng);
    ASSERT_EQ(sol, 0);
}

TEST(GainBucket, KLFMpass) {
    auto H =
    pmondriaan::read_hypergraph("../test/data/matrices/dolphins/dolphins.mtx", "degree")
    .value();
    std::mt19937 rng(1);
    pmondriaan::options opts;
    opts.KLFM_max_passes = 1;
    opts.metric = pmondriaan::m::cut_net;
    auto sol = pmondriaan::initial_partitioning(H, 163, 163, opts, rng);
    ASSERT_LE(H.weight_part(0), 163);
    ASSERT_LE(H.weight_part(1), 163);
    ASSERT_GE(sol, 18);
}

} // namespace
} // namespace pmondriaan
