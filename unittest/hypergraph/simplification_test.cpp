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

std::string mtx_three_size_edge = R"(%%MatrixMarket matrix coordinate real symmetric
3 3 3
1 1 1.0
1 2 1.0
1 3 1.0)";

TEST(Simplify, DuplicateNets) {
    std::stringstream mtx_ss(mtx_three_size_edge);
    auto hypergraph = read_hypergraph_istream(mtx_ss, "one");
    auto h3 = hypergraph.value();
    std::vector<long> verts;
    verts.push_back(h3.nets()[0].vertices()[0]);
    verts.push_back(h3.nets()[0].vertices()[1]);
    verts.push_back(h3.nets()[0].vertices()[2]);
    h3.add_net(1, verts, h3.nets()[0].cost());
    ASSERT_EQ(h3.nets().size(), 2);
    simplify_duplicate_nets(h3);
    ASSERT_EQ(h3.nets().size(), 1);
    ASSERT_EQ(h3.nets()[0].cost(), 2);
}

TEST(Simplify, BreakTriples) {
    std::stringstream mtx_ss(mtx_three_size_edge);
    auto hypergraph = read_hypergraph_istream(mtx_ss, "one");
    auto h3 = hypergraph.value();
    auto origsize1 = h3.vertices()[0].deg();
    auto origsize2 = h3.vertices()[1].deg();
    auto origsize3 = h3.vertices()[2].deg();
    break_triples(h3);
    ASSERT_LT(origsize1, h3.vertices()[0].deg());
    ASSERT_LT(origsize2, h3.vertices()[1].deg());
    ASSERT_LT(origsize3, h3.vertices()[2].deg());
}

} // namespace
} // namespace pmondriaan


