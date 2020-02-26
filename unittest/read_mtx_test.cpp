#include "pmondriaan.hpp"

#include "gtest/gtest.h"

namespace pmondriaan {

TEST(ReadMTX, ReadingFile) {
    auto H_missing = read_hypergraph("file_does_not_exist.mtx", "one");
    ASSERT_TRUE(!H_missing);

    auto H_exists = read_hypergraph("../test/data/matrices/cage3/cage3.mtx", "one");
    ASSERT_TRUE(H_exists);
}

} // namespace pmondriaan

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
