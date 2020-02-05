#include <iostream>
#include <fstream>
#include <vector>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include <pmondriaan.hpp>


int main () {	
	
	bulk::thread::environment env;

    env.spawn(env.available_processors(), [](bulk::world& world) {
		srand(world.rank() + 1);
		auto H = pmondriaan::read_hypergraph("../test/data/matrices/dolphins/dolphins.mtx", world, "degree");
		auto opts = pmondriaan::options();
		opts.sample_size = 5;
		opts.coarsening_max_clustersize = 5;
		opts.lp_max_iterations = 4;
		H.print();
		pmondriaan::coarsen_hypergraph(world, H, opts, "label propagation");
	});
	return 0;
}