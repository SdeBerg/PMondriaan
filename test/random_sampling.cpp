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
		auto hypergraph = pmondriaan::read_hypergraph("../test/data/matrices/gemat11/gemat11.mtx", world, "degree");
	
		pmondriaan::coarsen_hypergraph(world, hypergraph);
	});
	return 0;
}