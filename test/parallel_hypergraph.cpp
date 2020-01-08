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

int main() {
	
	environment env;
	
    env.spawn(env.available_processors(), [](bulk::world& world) {
        int s = world.rank();
        //int p = world.active_processors();

		auto hypergraph = pmondriaan::read_hypergraph("../test/data/matrices/cage3/cage3.mtx", world);
		
		for (auto i = 0u; i < hypergraph.size(); i++) {
			auto v = hypergraph(i);
			world.log("processor: %d, id: %d, weight: %d, \n nets: ", s, v.id(), v.weight());
			auto nets = v.nets();
			for (auto n : nets) {
				world.log("%d ", n);
			}
		}
    });

    return 0;
}
