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
		int s = world.rank();
		
		srand(world.rank() + 1);
		auto hypergraph = pmondriaan::read_hypergraph("../test/data/matrices/gemat11/gemat11.mtx", world, "one");
		if(!hypergraph) {
			std::cerr << "Error: failed to load hypergraph\n";
			return;
		}
		auto H = hypergraph.value();

		/*int count = 0;
		while (count < p) {
			if (s == count) { H.print(); }
			world.sync();
			count ++;
		}*/
		
		auto opts = pmondriaan::options();
		opts.sample_size = 200;
		opts.coarsening_max_clustersize = 400;
		opts.lp_max_iterations = 4;
		opts.coarsening_maxrounds = 15;
		opts.coarsening_nrvertices = 200;
		
		recursive_bisect(world, H, "multilevel", "label propagation", "cutnet", 2, 0.03, 0.03, opts);
		
		auto lb = pmondriaan::load_balance(world, H, 2);
		auto cutsize = pmondriaan::cutsize(world, H, "lambda1");
		
		/*for (auto& net : H.nets()) {
			if (s == 0) {
				world.log("net: %d", net.id());
			}
			int count = 0;
			while (count < p) {
				if (s == count) { 
					for (auto v : net.vertices()) {
						world.log("%d, part: %d", v, H(H.local_id(v)).part());
					}
				}
				world.sync();
				count ++;
			}
		}*/
		
		if (s == 0) {
			world.log("Load balance of partitioning found: %lf", lb);
			world.log("Cutsize of partitioning found: %d", cutsize);
		}
		world.sync();
	});
	return 0;
}