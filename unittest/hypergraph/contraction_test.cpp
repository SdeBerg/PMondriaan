#include "pmondriaan.hpp"

#include <random>

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

std::string test_mtx = R"(%%MatrixMarket matrix coordinate real general
3 4 5
1 1 1.0
2 1 1.0
2 2 1.0
1 3 1.0
3 4 1.0
)";

TEST(Contraction, MergeFreeVertices) {
    std::stringstream mtx_ss(test_mtx);
    auto hypergraph = read_hypergraph_istream(mtx_ss, "one");
    auto H = hypergraph.value();
    remove_free_nets(H, 1);
    ASSERT_EQ(H.size(), 4);
    ASSERT_EQ(H.global_size(), 4);
    auto C = pmondriaan::contraction();
    C.merge_free_vertices(H);
    ASSERT_EQ(H.size(), 3);
    ASSERT_EQ(H.global_size(), 3);
    H(0).set_part(0);
    H(1).set_part(1);
    H(2).set_part(0);
    std::mt19937 rng(1);
    auto result = C.assign_free_vertices(H, 3, 3, rng);
    ASSERT_EQ(result[0], 2);
    ASSERT_EQ(result[1], 2);
}

TEST(Contraction, ParallelMergeFreeVertices) {
    environment env;
    env.spawn(2, [](bulk::world& world) {
        std::stringstream mtx_ss(test_mtx);
        auto hypergraph = read_hypergraph_istream(mtx_ss, world, "one");
        auto H = hypergraph.value();
        remove_free_nets(world, H, 1);
        ASSERT_EQ(H.global_size(), 4);
        auto C = pmondriaan::contraction();
        C.merge_free_vertices(world, H);
        ASSERT_EQ(H.global_size(), 3);
        if (world.rank() == 0) {
            H(0).set_part(0);
            H(1).set_part(1);
        }
        if (world.rank() == 1) {
            H(0).set_part(0);
        }
        std::mt19937 rng(1);
        auto result = C.assign_free_vertices(world, H, 3, 3, rng);
        ASSERT_EQ(result[0], 2);
        ASSERT_EQ(result[1], 2);
    });
}

} // namespace
} // namespace pmondriaan
