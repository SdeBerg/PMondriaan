#include "pmondriaan.hpp"

#include "gtest/gtest.h"
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
using environment = bulk::mpi::environment;
#else
#include <bulk/backends/thread/thread.hpp>
using environment = bulk::thread::environment;
#endif

namespace pmondriaan {

std::string mtx_three_nonzeros = R"(%%MatrixMarket matrix coordinate real general
3 3 4
1 1 1.0
2 1 1.0
2 2 1.0
1 3 1.0
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

TEST(Cutsize, CutsizeCutnet) {
    environment env;
    env.spawn(2, [](bulk::world& world) {
        std::stringstream mtx_ss(mtx_three_nonzeros);
        auto hypergraph = read_hypergraph_istream(mtx_ss, world, "degree");
        auto H = hypergraph.value();
        world.log("H(0) %d", H(0).id());
        world.log("global1 %d", H.global_number_nets());
        world.log("global size %d", H.global_size());
        world.sync();
        /*if (world.rank() == 0) {
            H(0).set_part(0);
            H(1).set_part(1);
        }
        if (world.rank() == 1) {
            H(0).set_part(1);
        }
        auto cut = pmondriaan::cutsize(world, H, pmondriaan::m::cut_net);
        if (world.rank() == 0) {
            ASSERT_EQ(cut, 4);
        }
        */
    });
}

TEST(Bisect, BisectLP) {
    auto H =
    pmondriaan::read_hypergraph("../test/data/matrices/dolphins/dolphins.mtx", "one")
    .value();
    std::mt19937 rng(1);
    auto C = std::vector<std::vector<long>>(H.nets().size(), std::vector<long>(2, 0));
    auto L = label_propagation_bisect(H, C, 100, 40, 40, rng);
    ASSERT_EQ(L.size(), H.size());
}

TEST(GainBucket, GainBucketInit) {
    std::stringstream mtx_ss(mtx_three_nonzeros);
    auto H = read_hypergraph_istream(mtx_ss, "one").value();
    H(0).set_part(0);
    H(1).set_part(1);
    H(2).set_part(0);
    auto C = std::vector<std::vector<long>>(3);
    C[0] = {2, 0};
    C[1] = {1, 1};
    C[2] = {0, 0};
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
    pmondriaan::read_hypergraph("../test/data/matrices/dolphins/dolphins.mtx", "one")
    .value();
    std::mt19937 rng(1);
    pmondriaan::options opts;
    opts.KLFM_max_passes = 1;
    opts.metric = pmondriaan::m::cut_net;
    auto sol = pmondriaan::initial_partitioning(H, 35, 35, opts, rng);
    ASSERT_LE(H.weight_part(0), 35);
    ASSERT_LE(H.weight_part(1), 35);
    ASSERT_GE(sol, 6);
}

} // namespace pmondriaan

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
