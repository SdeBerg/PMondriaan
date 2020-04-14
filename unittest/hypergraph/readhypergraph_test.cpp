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

TEST(ReadHypergraph, read_hypergraph) {
    auto H_missing = read_hypergraph("file_does_not_exist.mtx", "one");
    ASSERT_TRUE(!H_missing);

    auto H_exists =
    read_hypergraph("../test/data/matrices/cage3/cage3.mtx", "one");
    ASSERT_TRUE(H_exists);

    std::stringstream mtx_ss(mtx_three_nonzeros);
    auto H = read_hypergraph_istream(mtx_ss, "degree");
    ASSERT_TRUE(H);
    ASSERT_EQ(H->size(), 3);
}

} // namespace
} // namespace pmondriaan


