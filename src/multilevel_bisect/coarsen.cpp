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
	//auto indices_samples = sample_random(HC, options.sample_size());
	auto indices_samples = sample_random(H, 1);
	
	//we now send the samples and the processor id to all processors
	auto sample_queue = bulk::queue<int, int, int[]>(world);
	for (auto i = 0u; i < indices_samples.size(); i++) {
		for (int t = 0; t < p; t++) {
			sample_queue(t).send(s, H(indices_samples[i]).id(), H(indices_samples[i]).nets());
		}
	}
	world.sync();
	
	auto number_samples = std::vector<int>(p, 0);
	//auto ip = std::vector<std::vector<int>>(H.size(), std::vector<int>(p * options.sample_size(), 0));
	auto ip = std::vector<std::vector<int>>(H.size(), std::vector<int>(p * 1, 0));
	
	//we store the ids of the received samples
	auto ids = std::vector<std::vector<int>>(p, std::vector<int>(options.sample_size(), 0));
	
	for (auto& [t, sample_id, sample_nets] : std::move(sample_queue)) {
		ids[t][number_samples[t]] = sample_id;
		for (auto n_id : sample_nets) {
			for (auto u_id : H.net(n_id).vertices()) {
				//ip[H.local_id(u_id)][t * options.sample_size + samples[t]]++;
				ip[H.local_id(u_id)][t * 1 + number_samples[t]]++;
			}
		}
		number_samples[t]++;
	}
	
	
	

	return H;
}

} // namespace pmondriaan
