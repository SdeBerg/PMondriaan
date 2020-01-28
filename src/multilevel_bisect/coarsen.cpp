#include <vector>
#include <array>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "bisect.hpp"
#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/sample.hpp"

namespace pmondriaan {

/**
 * Coarses the hypergraph H and returns a hypergraph HC.
 */
pmondriaan::hypergraph coarsen_hypergraph(bulk::world& world, pmondriaan::hypergraph& H) {
	
	int s = world.rank();
	int p = world.active_processors();
	
	//we first select ns samples	
	//auto indices_samples = sample_random(HC, options.sample_size);
	auto indices_samples = sample_random(H, 10);
	
	//we now send the samples and the processor id to all processors
	auto sample_queue = bulk::queue<int, int, int[]>(world);
	for (auto i = 0u; i < indices_samples.size(); i++) {
		for (int t = 0; t < p; t++) {
			sample_queue(t).send(s, indices_samples[i], H(indices_samples[i]).nets());
		}
	}
	world.sync();
	
	/*auto samples = std::vector<std::vector<int[]>>(p);
	auto ip = std::array<std::array<int, p * options.sample_size>, H.size()>();
	
		for (const auto& [id, index, sample] : std::move(sample_queue)) {
			world.log("t: %d, index: %d, sample[0]: %d", id, index, sample[0]);
			samples[id].push_back(sample);
			for (auto n : sample) {
				for (auto u : H.net(n).vertices()) {
					ip[H.global_id(u.id())][nr_sample]++;
				}
			}
		}*/
	

	return H;
}

} // namespace pmondriaan
