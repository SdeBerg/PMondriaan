#include <iostream>

#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
using environment = bulk::mpi::environment;
#else
#include <bulk/backends/thread/thread.hpp>
using environment = bulk::thread::environment;
#endif

#include <pmondriaan.hpp>

int main() {
	
	environment env;
	
    env.spawn(env.available_processors(), [](bulk::world& world) {
        int s = world.rank();
        //int p = world.active_processors();
		
		srand(s + 1);

		auto hypergraph = pmondriaan::read_hypergraph("../test/data/matrices/west0381/west0381.mtx", world);
		
		auto weights = bisect_random(world, hypergraph, 0.05);
		world.log("weight part 0: %d, weight part 1: %d", weights[0], weights[1]);
		auto eps = compute_load_balance(world, hypergraph, 2);
		world.log("load imbalance: %lf", eps);
	
	});

    return 0;
}
