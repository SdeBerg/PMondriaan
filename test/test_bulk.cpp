#include <pmondriaan.hpp>

int main() {
    bulk::thread::environment env;

    env.spawn(env.available_processors(), [](bulk::world& world) {
        int s = world.rank();
        int p = world.active_processors();
		
		auto vertices_queue = bulk::queue<int, int[]>(world);

        auto hypergraph = pmondriaan::read_hypergraph("../test/data/matrices/west0381/west0381.mtx", world);
		
		if(s == 0) {
			hypergraph(0).set_id(10);
			for(int t = 0; t < p; t++) {
				vertices_queue(t).send(hypergraph(0).id(), hypergraph(0).nets());
			}
			hypergraph(0).set_id(5);
			hypergraph.vertices().erase(hypergraph.vertices().begin() + 1);
		}
		
		world.sync();
		
		for (const auto& [index, nets] : vertices_queue) {
			world.log("%d", index);
			world.log("%d", nets[1]);
		}
		
    });

    return 0;
}
