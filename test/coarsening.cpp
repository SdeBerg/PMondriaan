#include <fstream>
#include <iostream>
#include <vector>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include <pmondriaan.hpp>


int main() {

    bulk::thread::environment env;

    env.spawn(env.available_processors(), [](bulk::world& world) {
        int p = world.active_processors();
        int s = world.rank();

        srand(world.rank() + 1);
        auto H = pmondriaan::read_hypergraph(
        "../test/data/matrices/dolphins/dolphins.mtx", world, "degree");

        auto opts = pmondriaan::options();
        opts.sample_size = 5;
        opts.coarsening_max_clustersize = 20;
        opts.lp_max_iterations = 4;
        opts.coarsening_maxrounds = 10;

        auto HC_new = pmondriaan::create_new_hypergraph(world, H, 0, H.size());
        auto C = pmondriaan::contraction();
        auto HC = pmondriaan::coarsen_hypergraph(world, HC_new, C, opts, "label propagation");

        for (auto& net : HC.nets()) {
            if (s == 0) {
                world.log("net: %d", net.id());
            }
            int count = 0;
            while (count < p) {
                if (s == count) {
                    for (auto v : net.vertices()) {
                        world.log("%d", v);
                    }
                }
                world.sync();
                count++;
            }
        }
    });
    return 0;
}