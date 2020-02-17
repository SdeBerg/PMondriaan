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
		//int p = world.active_processors();
		//int s = world.rank();
		
		srand(world.rank() + 1);
		auto H = pmondriaan::read_hypergraph("../test/data/matrices/west0381/west0381.mtx", world, "degree");

		/*int count = 0;
		while (count < p) {
			if (s == count) { H.print(); }
			world.sync();
			count ++;
		}*/
		
		auto opts = pmondriaan::options();
		opts.sample_size = 10;
		opts.coarsening_max_clustersize = 30;
		opts.lp_max_iterations = 4;
		opts.coarsening_maxrounds = 10;

		//recursive_bisect(world, H, "multilevel", "label propagation", "cutnet", 2, 0.03, 0.03, opts);
		auto HC_new = pmondriaan::create_new_hypergraph(world, H, 0, H.size());
		auto HC = pmondriaan::coarsen_hypergraph(world, HC_new, opts, "label propagation");
		world.log("global size: %d", HC.global_size()); world.sync();
		/*count = 0;
		while (count < p) {
			if (s == count) { HC.print(); }
			world.sync();
			count ++;
		}
		
		for (auto net : HC.nets()) {
			world.log("net: %d", net.id());
			for (auto v : net.vertices()) {
				world.log("%d", v);
			}
			world.sync();
		}*/
	});
	return 0;
}