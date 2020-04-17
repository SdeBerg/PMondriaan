#include "pmondriaan.hpp"

#include "gtest/gtest.h"

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
using environment = bulk::mpi::environment;
#else
#include <bulk/backends/thread/thread.hpp>
using environment = bulk::thread::environment;
#endif

namespace pmondriaan {
namespace {

std::string mtx_three_nonzeros = R"(%%MatrixMarket matrix coordinate real general
3 3 4
1 1 1.0
2 1 1.0
2 2 1.0
1 3 1.0
)";

TEST(KLFMParallel, InitCounts) {
    environment env;
    env.spawn(2, [](bulk::world& world) {
        std::stringstream mtx_ss(mtx_three_nonzeros);
        auto hypergraph = read_hypergraph_istream(mtx_ss, world, "degree");
        auto H = hypergraph.value();
        if (world.rank() == 0) {
            H(0).set_part(0);
            H(1).set_part(1);
        }
        if (world.rank() == 1) {
            H(0).set_part(1);
        }
        auto C = pmondriaan::init_counts(world, H);
        ASSERT_EQ(C[0][0], 1);
        ASSERT_EQ(C[0][1], 1);
        if (world.rank() == 0) {
            ASSERT_EQ(C[1][0], 1);
            ASSERT_EQ(C[1][1], 1);
        }
    });
}

TEST(KLFMParallel, InitPreviousC) {
    environment env;
    env.spawn(2, [](bulk::world& world) {
        auto s = world.rank();
        std::stringstream mtx_ss(mtx_three_nonzeros);
        auto hypergraph = read_hypergraph_istream(mtx_ss, world, "degree");
        auto H = hypergraph.value();
        if (s == 0) {
            H(0).set_part(0);
            H(1).set_part(1);
        }
        if (s == 1) {
            H(0).set_part(1);
        }
        auto C = pmondriaan::init_counts(world, H);
        ASSERT_EQ(C[0][0], 1);
        ASSERT_EQ(C[0][1], 1);
        if (s == 0) {
            ASSERT_EQ(C[1][0], 1);
            ASSERT_EQ(C[1][1], 1);
        }
        auto net_partition = bulk::block_partitioning<1>({H.global_number_nets()}, 2);
        auto previous_C = bulk::coarray<long>(world, net_partition.local_count(s) * 2);
        auto cost_my_nets = bulk::coarray<long>(world, net_partition.local_count(s));
        auto cut_size_my_nets =
        init_previous_C(world, H, C, previous_C, cost_my_nets, net_partition);
        if (s == 0) {
            ASSERT_EQ(cut_size_my_nets, 2);
        }
        if (s == 1) {
            ASSERT_EQ(cut_size_my_nets, 0);
        }
    });
}

} // namespace
} // namespace pmondriaan
