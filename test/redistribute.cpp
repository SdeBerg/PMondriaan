#include <pmondriaan.hpp>

int main() {
    bulk::thread::environment env;

    env.spawn(env.available_processors(), [](bulk::world& world) {
        //int s = world.rank();
        //int p = world.active_processors();
		
		auto H = pmondriaan::read_hypergraph("../test/data/matrices/west0381/west0381.mtx", world, "degree");

		recursive_bisect(world, H, "random", 4, 0.05, 0.05);
		/*for (auto& v : H.vertices()) {
			if (s == 0 || s == 1) {
				if (v.part() == 1) {
					world.log("wrong");
				}
			}
			else {
				if (v.part() == 0) {
					world.log("wrong");
				}
			}
		}*/
		
    });

    return 0;
}
