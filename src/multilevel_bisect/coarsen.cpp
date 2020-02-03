#include <vector>
#include <array>
#include <algorithm>

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
	
	/*int s = world.rank();
	int p = world.active_processors();
	
	//we first select ns samples	
	auto indices_samples = sample_random(HC, options.sample_size());
	//auto indices_samples = sample_random(H, 1);
	
	//we now send the samples and the processor id to all processors
	auto sample_queue = bulk::queue<int, int, int[]>(world);
	for (auto i = 0u; i < indices_samples.size(); i++) {
		for (int t = 0; t < p; t++) {
			sample_queue(t).send(s, i, H(indices_samples[i]).nets());
		}
	}
	world.sync();
	
	int total_samples = p * options.sample_size();
	
	auto ip = std::vector<std::vector<double>>(H.size(), std::vector<double>(total_samples, 0.0));
	
	auto degree_samples = std::vector<int>(total_samples);
	
	for (auto& [t, number_sample, sample_nets] : std::move(sample_queue)) {
		degree_samples[t * options.sample_size() + number_sample] = sample_nets.size();
		for (auto n_id : sample_nets) {
			for (auto u_id : H.net(n_id).vertices()) {
				ip[H.local_id(u_id)][t * options.sample_size() + number_sample] += 1.0/(double)H.net(n_id).size();
				//ip[H.local_id(u_id)][t * 1 + number_sample]+= 1.0/(double)H.net(n_id).size();
			}
		}
	}
	
	//we set the ip of all local samples with itself to 0, so they will not match themselves
	for (auto i = 0u; i < indices_samples.size(); i++) {
		ip[indices_samples[i]][s * options.sample_size() + i] = 0;
	}
	
	//find best sample for vertex v and add it to the list of that sample
	auto requested_matches = std::vector<std::vector<std::pair<int, double>>>(total_samples);
	for (auto& v : H.vertices()) {
		double max_ip = 0.0;
		int best_match;
		auto local_id = H.local_id(v.id());
		for (auto u = 0; u < total_samples; u++) {
			double ip_vu = ip[local_id][u];
			if (ip_vu > 0) {
				ip_vu *= 1.0/(double)std::min(v.degree(), degree_samples[u]);
				if (ip_vu > max_ip) {
					max_ip = ip_vu;
					best_match = u;
				}
			}
		}
		requested_matches[best_match].push_back(std::make_pair(local_id, max_ip));
	}
	
	for (auto& match_list : requested_matches) {
		std::sort(match_list.begin(), match_list.end(), [] (std::pair const& match1, std::pair const& ,match2) { return match1.second > match2.second; });
	}*/
	
	
	

	return H;
}

} // namespace pmondriaan
