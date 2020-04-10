#include "pmondriaan.hpp"

#include "gtest/gtest.h"

namespace pmondriaan {
  namespace {

TEST(Bisect, BisectLP) {
    auto H =
    pmondriaan::read_hypergraph("../test/data/matrices/dolphins/dolphins.mtx", "one")
    .value();
    std::mt19937 rng(1);
    auto C = std::vector<std::vector<long>>(H.nets().size(), std::vector<long>(2, 0));
    auto L = label_propagation_bisect(H, C, 100, 40, 40, rng);
    ASSERT_EQ(L.size(), H.size());
}

} // namespace
} // namespace mondriaan
