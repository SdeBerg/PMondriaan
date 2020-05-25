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

TEST(KLFMParallel, KLFMParPass) {
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

        pmondriaan::options opts;
        opts.KLFM_max_passes = 1;
        opts.KLFM_par_number_send_moves = 1;
        opts.metric = pmondriaan::m::cut_net;
        std::mt19937 rng(world.rank() + 1);
        auto cut = KLFM_par(world, H, C, global_weight_part(world, H, 0),
                            global_weight_part(world, H, 1), 3, 3, opts, rng);
        ASSERT_LE(cut, 2);
        ASSERT_EQ(cut, pmondriaan::cutsize(world, H, opts.metric));
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

TEST(KLFMParallel, RejectMoves) {
    environment env;
    env.spawn(2, [](bulk::world& world) {
        auto s = world.rank();
        // Stores the proposed moves as: gain, weight change, processor id
        auto moves_queue = bulk::queue<long, long, int, long>(world);
        if (s == 0) {
            moves_queue(0).send(5, -3, s, 0);
            moves_queue(0).send(3, -4, s, 1);
            moves_queue(0).send(1, 3, s, 2);
            moves_queue(0).send(-2, -1, s, 3);
            moves_queue(0).send(-2, 2, s, 4);
        } else {
            moves_queue(0).send(5, -1, s, 0);
            moves_queue(0).send(4, -2, s, 1);
            moves_queue(0).send(2, -3, s, 2);
            moves_queue(0).send(1, 4, s, 3);
            moves_queue(0).send(-5, 0, s, 4);
        }
        world.sync();
        auto prev_total_weights = std::array<long, 2>({20, 18});
        long max_weight_0 = 22;
        long max_weight_1 = 22;
        auto rejected = reject_unbalanced_moves(world, moves_queue, prev_total_weights,
                                                max_weight_0, max_weight_1);
        if (s == 0) {
            ASSERT_EQ(rejected, -1);
        }
        if (s == 1) {
            ASSERT_EQ(rejected, 0);
        }
    });
}

TEST(KLFMParallel, KLFMPar) {
    environment env;
    env.spawn(2, [](bulk::world& world) {
        auto hypergraph =
        read_hypergraph("../test/data/matrices/dolphins/dolphins.mtx", world, "degree");
        auto H = hypergraph.value();
        pmondriaan::interval labels = {0, 1};
        std::mt19937 rng(world.rank() + 1);
        bisect_random(H, 85, 85, 0, H.size(), labels, rng);
        ASSERT_LE(global_weight_part(world, H, 0), 170);
        ASSERT_LE(global_weight_part(world, H, 1), 170);

        pmondriaan::options opts;
        opts.KLFM_max_passes = 1;
        opts.KLFM_par_number_send_moves = 3;
        opts.metric = pmondriaan::m::cut_net;
        auto C = init_counts(world, H);
        auto cut = KLFM_par(world, H, C, global_weight_part(world, H, 0),
                            global_weight_part(world, H, 1), 170, 170, opts, rng);
        ASSERT_EQ(cut, pmondriaan::cutsize(world, H, opts.metric));
        ASSERT_LE(global_weight_part(world, H, 0), 170);
        ASSERT_LE(global_weight_part(world, H, 1), 170);
    });
} // namespace

} // namespace
} // namespace pmondriaan
