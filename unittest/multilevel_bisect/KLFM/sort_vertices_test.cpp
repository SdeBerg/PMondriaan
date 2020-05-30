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

TEST(KLFMSortVertices, KLFMSortVerticesInit) {
    std::stringstream mtx_ss(mtx_three_nonzeros);
    auto hypergraph = read_hypergraph_istream(mtx_ss, "degree");
    auto H = hypergraph.value();
    H(0).set_part(1);
    H(1).set_part(0);
    H(2).set_part(1);
    auto C = pmondriaan::init_counts(H);
    ASSERT_EQ(C[0][0], 0);
    ASSERT_EQ(C[0][1], 2);
    ASSERT_EQ(C[1][0], 1);
    ASSERT_EQ(C[1][1], 1);

    H.sort_vertices_on_part(C);
    ASSERT_EQ(H.net(1).vertices()[0], 1);
    ASSERT_EQ(H.net(1).vertices()[1], 0);

    auto H2 =
    read_hypergraph("../test/data/matrices/dolphins/dolphins.mtx", "degree").value();
    pmondriaan::interval labels = {0, 1};
    std::mt19937 rng(1);
    bisect_random(H2, 163, 163, 0, H2.size(), labels, rng);
    auto C2 = pmondriaan::init_counts(H2);

    H2.sort_vertices_on_part(C2);
    for (auto net : H2.nets()) {
        size_t index = 0;
        while (index < net.size() && H2(H2.local_id(net.vertices()[index])).part() == 0) {
            index++;
        }
        while (index < net.size()) {
            ASSERT_EQ(H2(H2.local_id(net.vertices()[index])).part(), 1);
            index++;
        }
    }

    H2.move(10, C2);
    H2.move(13, C2);
    H2.move(30, C2);
    H2.move(55, C2);

    for (auto net : H2.nets()) {
        size_t index = 0;
        while (index < net.size() && H2(H2.local_id(net.vertices()[index])).part() == 0) {
            index++;
        }
        while (index < net.size()) {
            ASSERT_EQ(H2(H2.local_id(net.vertices()[index])).part(), 1);
            index++;
        }
    }
}

} // namespace
} // namespace pmondriaan
