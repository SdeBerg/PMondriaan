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
		opts.sample_size = 5;
		opts.coarsening_max_clustersize = 20;
		opts.lp_max_iterations = 4;
		opts.coarsening_maxrounds = 10;

		auto HC_new = pmondriaan::create_new_hypergraph(world, H, 0, H.size());
		auto C = pmondriaan::contraction();
		auto HC = pmondriaan::coarsen_hypergraph(world, HC_new, C, opts, "label propagation");
		world.log("global size: %d and weight: %d", H.global_size(), H.global_weight(world)); world.sync();
		world.log("global size: %d and weight: %d", HC.global_size(), HC.global_weight(world)); world.sync();
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