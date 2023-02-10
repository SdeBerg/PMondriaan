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
5 4 7
1 1 1.0
1 3 1.0
2 3 1.0
2 1 1.0
3 4 1.0
4 2 1.0
5 2 1.0
)";

TEST(Simplify, SimplifyDuplicateNets) {
    std::stringstream mtx_ss(test_mtx);
    auto hypergraph = read_hypergraph_istream(mtx_ss, "one");
    auto H = hypergraph.value();
    ASSERT_EQ(H.nets().size(), 5);
    simplify_duplicate_nets(H);
    ASSERT_EQ(H.nets().size(), 3);
}

} // namespace
} // namespace pmondriaan
