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

TEST(WritePartitioning, StartParts) {
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
            H(0).set_part(0);
        }
        auto result = pmondriaan::start_parts(world, H, 2);
        ASSERT_EQ(result[0], 0);
        ASSERT_EQ(result[1], 3);
        ASSERT_EQ(result[2], 4);
    });
}

} // namespace
} // namespace pmondriaan
