#include "pmondriaan.hpp"

#include "gtest/gtest.h"

namespace pmondriaan {

std::string mtx_three_nonzeros = R"(%%MatrixMarket matrix coordinate real general
3 3 3
1 1 1.0
2 1 1.0
3 1 1.0
)";

TEST(ReadMTX, ReadingFile) {
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

TEST(BisectLP, LPBisecting) {
    auto H =
    pmondriaan::read_hypergraph("../test/data/matrices/dolphins/dolphins.mtx", "one")
    .value();
    std::mt19937 rng(1);
    auto L = label_propagation_bisect(H, 100, 40, 40, rng);
    ASSERT_EQ(L.size(), H.size());
}

} // namespace pmondriaan

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
