#include <iostream>
#include <fstream>


#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
using environment = bulk::mpi::environment;
#else
#include <bulk/backends/thread/thread.hpp>
using environment = bulk::thread::environment;
#endif

int main() {
	
	environment env;
	
    env.spawn(4, [](bulk::world& world) {
        int s = world.rank();
        //int p = world.active_processors();

		auto partitioning = bulk::block_partitioning<1>({5}, {4});
		world.log("s: %d, local_count: %d", s, partitioning.local_count(s));
		
		world.log("v: %d, owner: %d", 3, partitioning.owner({3}));
		world.log("s: %d, v_glob: %d, v_loc: %d", s, 3, partitioning.local({3})[0]);
		world.log("s: %d, v_loc: %d, v_glob: %d", s, 0, partitioning.global({0}, s)[0]);
    });

    return 0;
}
