#include <array>
#include <fstream>
#include <iostream>
#include <random>
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

    env.spawn(4, [](bulk::world& world) {
        // int p = world.active_processors();
        int s = world.rank();

        srand(world.rank() + 1);
        auto hypergraph = pmondriaan::read_hypergraph(
        "../test/data/matrices/dolphins/dolphins.mtx", world, "one");
        if (!hypergraph) {
            std::cerr << "Error: failed to load hypergraph\n";
            return;
        }
        auto H = hypergraph.value();

        auto opts = pmondriaan::options();
        opts.sample_size = 200;
        opts.coarsening_max_clustersize = 400;
        opts.lp_max_iterations = 4;
        opts.coarsening_maxrounds = 15;
        opts.coarsening_nrvertices = 200;
        opts.KLFM_max_passes = 1;
        opts.KLFM_par_number_send_moves = 3;

        std::mt19937 rng(s + 1);
        pmondriaan::interval i = {0, 1};
        bisect_random(H, 18, 18, 0, H.size(), i, rng);
        auto counts = pmondriaan::init_counts(H);
        /*auto weights = std::array<long, 2>();
        weights[0] = H.weight_part(0);
        weights[1] = H.weight_part(1);
        KLFM_pass(H, counts, pmondriaan::cutsize(H, pmondriaan::m::cut_net),
                  weights, 35, 35, rng);*/
        KLFM_par(world, H, counts, H.weight_part(0), H.weight_part(1), 35, 35, opts, rng);


        /*recursive_bisect(world, H, "multilevel", "label propagation",
        "cutnet", 2, 0.03, 0.03, opts);

        auto lb = pmondriaan::load_balance(world, H, 2);
        auto cutsize = pmondriaan::cutsize(world, H, "lambda1");

        if (s == 0) {
            world.log("Load balance of partitioning found: %lf", lb);
            world.log("Cutsize of partitioning found: %d", cutsize);
        }
        world.sync();*/
    });
    return 0;
}