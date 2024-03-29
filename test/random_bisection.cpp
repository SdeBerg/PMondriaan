#include <iostream>
#include <fstream>

#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
using environment = bulk::mpi::environment;
#else
#include <bulk/backends/thread/thread.hpp>
using environment = bulk::thread::environment;
#endif

#include <hypergraph/readhypergraph.hpp>
#include <bisect.hpp>

int main() {
	
	environment env;
	
    env.spawn(env.available_processors(), [](bulk::world& world) {
        int s = world.rank();
        //int p = world.active_processors();
		
		srand(s + 1);

		auto hypergraph = pmondriaan::read_hypergraph("../test/data/matrices/cage3/cage3.mtx", world);
		
		auto sizes = pmondriaan::global_net_sizes(world, hypergraph);
		
		if (s == 0) {
			for (auto s : sizes) {
				world.log("%d", s);
			}
		}
		
		auto weights = bisect_random(world, hypergraph, 0.05);
		world.log("weight part 0: %d, weight part 1: %d", weights[0], weights[1]);
		
    });

    return 0;
}
