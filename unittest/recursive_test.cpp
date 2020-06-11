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

TEST(RecursiveBisect, SeqRecursiveBisect) {
    environment env;
    env.spawn(1, [](bulk::world& world) {
        auto H = pmondriaan::read_hypergraph(
                 "../test/data/matrices/dolphins/dolphins.mtx", "degree")
                 .value();
        auto old_size = H.size();
        auto old_nets = H.nets().size();
        std::mt19937 rng(1);
        pmondriaan::options opts;
        opts.KLFM_max_passes = 10;
        opts.metric = pmondriaan::m::cut_net;
        opts.bisection_mode = pmondriaan::bisection::multilevel;
        opts.coarsening_max_clustersize = 5;
        opts.lp_max_iterations = 10;
        opts.coarsening_nrvertices = 40;
        opts.coarsening_maxrounds = 1;
        recursive_bisect(world, H, 7, 0.1, 0.1, opts);
        auto lb = pmondriaan::load_balance(world, H, 7);
        ASSERT_LE(lb, 0.1);
        for (auto v : H.vertices()) {
            ASSERT_NE(v.part(), -1);
        }
        ASSERT_EQ(old_size, H.size());
        ASSERT_EQ(old_nets, H.nets().size());
    });
}

TEST(RecursiveBisect, ParRecursiveBisect) {
    environment env;
    env.spawn(3, [](bulk::world& world) {
        auto H = pmondriaan::read_hypergraph(
                 "../test/data/matrices/dolphins/dolphins.mtx", world, "degree")
                 .value();
        auto old_size = bulk::sum(world, H.size());
        auto local_old_nets = bulk::coarray<int>(world, H.global_number_nets());
        for (auto i = 0u; i < H.global_number_nets(); i++) {
            local_old_nets[i] = 0;
        }
        for (auto net : H.nets()) {
            local_old_nets[net.id()] = net.size();
        }
        auto old_nets =
        bulk::foldl_each(local_old_nets, [](auto& lhs, auto rhs) { lhs += rhs; });
        pmondriaan::options opts;
        opts.sample_size = 300;
        opts.KLFM_max_passes = 10;
        opts.metric = pmondriaan::m::cut_net;
        opts.sampling_mode = pmondriaan::sampling::label_propagation;
        opts.bisection_mode = pmondriaan::bisection::multilevel;
        opts.coarsening_max_clustersize = 5;
        opts.lp_max_iterations = 10;
        opts.coarsening_nrvertices = 40;
        opts.coarsening_maxrounds = 1;
        opts.KLFM_par_number_send_moves = 4;
        recursive_bisect(world, H, 7, 0.1, 0.1, opts);
        auto lb = pmondriaan::load_balance(world, H, 7);
        auto new_size = bulk::sum(world, H.size());
        ASSERT_EQ(old_size, new_size);
        for (auto v : H.vertices()) {
            ASSERT_NE(v.part(), -1);
        }
        auto local_new_nets = bulk::coarray<int>(world, H.global_number_nets());
        for (auto i = 0u; i < H.global_number_nets(); i++) {
            local_new_nets[i] = 0;
        }
        for (auto net : H.nets()) {
            local_new_nets[net.id()] = net.size();
        }
        auto new_nets =
        bulk::foldl_each(local_new_nets, [](auto& lhs, auto rhs) { lhs += rhs; });
        ASSERT_EQ(old_nets.size(), new_nets.size());
        for (auto i = 0u; i < old_nets.size(); i++) {
            ASSERT_EQ(old_nets[i], new_nets[i]);
        }
        ASSERT_LE(lb, 0.1);
    });
}

} // namespace
} // namespace pmondriaan
