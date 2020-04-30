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

TEST(BisectMultilevel, SeqBisectMultilevel) {
    environment env;
    env.spawn(1, [](bulk::world& world) {
        auto H = pmondriaan::read_hypergraph(
                 "../test/data/matrices/dolphins/dolphins.mtx", "degree")
                 .value();
        std::mt19937 rng(1);
        pmondriaan::options opts;
        opts.KLFM_max_passes = 10;
        opts.metric = pmondriaan::m::cut_net;
        opts.coarsening_max_clustersize = 5;
        opts.lp_max_iterations = 10;
        opts.coarsening_nrvertices = 40;
        opts.coarsening_maxrounds = 1;
        pmondriaan::interval labels = {0, 3};
        bisect_multilevel(world, H, opts, 163, 163, 0, H.size(), labels, rng);
        ASSERT_LE(H.weight_part(0), 163);
        ASSERT_LE(H.weight_part(3), 163);
        for (auto v : H.vertices()) {
            ASSERT_NE(v.part(), -1);
        }
        ASSERT_GE(pmondriaan::cutsize(H, opts.metric), 18);
        world.log("Cut %d", pmondriaan::cutsize(H, opts.metric));
    });
}

TEST(BisectMultilevel, ParBisectMultilevel) {
    environment env;
    env.spawn(3, [](bulk::world& world) {
        auto H = pmondriaan::read_hypergraph(
                 "../test/data/matrices/dolphins/dolphins.mtx", world, "degree")
                 .value();
        std::mt19937 rng(world.rank() + 1);
        pmondriaan::options opts;
        opts.sample_size = 4;
        opts.KLFM_max_passes = 10;
        opts.metric = pmondriaan::m::cut_net;
        opts.sampling_mode = pmondriaan::sampling::label_propagation;
        opts.coarsening_max_clustersize = 5;
        opts.lp_max_iterations = 10;
        opts.coarsening_nrvertices = 40;
        opts.coarsening_maxrounds = 1;
        opts.KLFM_par_number_send_moves = 4;
        pmondriaan::interval labels = {0, 3};
        bisect_multilevel(world, H, opts, 163, 163, 0, H.size(), labels, rng);
        ASSERT_LE(global_weight_part(world, H, 0), 163);
        ASSERT_LE(global_weight_part(world, H, 3), 163);
        ASSERT_EQ(global_weight_part(world, H, 1), 0);
        for (auto v : H.vertices()) {
            ASSERT_NE(v.part(), -1);
        }
        ASSERT_GE(pmondriaan::cutsize(world, H, opts.metric), 18);
        world.log("Cut %d", pmondriaan::cutsize(world, H, opts.metric));
    });
}

} // namespace
} // namespace pmondriaan
