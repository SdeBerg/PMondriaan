#include <pmondriaan.hpp>

int main() {
    bulk::thread::environment env;

    env.spawn(4, [](bulk::world& world) {
        //int s = world.rank();
        //int p = world.active_processors();
		
		auto H = pmondriaan::read_hypergraph("../test/data/matrices/west0381/west0381.mtx", world, "degree");
		
		recursive_bisect(world, H, "random", 4, 0.05, 0.05);
		
    });

    return 0;
}
