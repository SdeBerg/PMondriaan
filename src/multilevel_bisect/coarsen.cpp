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
pmondriaan::hypergraph coarsen_hypergraph(bulk::world& world, pmondriaan::hypergraph& H, pmondriaan::options& opts, std::string sampling_mode) {
	
	int s = world.rank();
	int p = world.active_processors();
	
	/* We first select ns samples */
	auto indices_samples = std::vector<int>();
	if (sampling_mode == "random") {
		indices_samples = sample_random(H, opts.sample_size);
	}
	else if (sampling_mode == "label propagation") {
		indices_samples = sample_lp(H, opts);
	}
	
	for (auto i : indices_samples) {
		world.log("s: %d, sample: %d", s, H(i).id());
	}
	world.sync();
	
	/* We now send the samples and the processor id to all processors */
	auto sample_queue = bulk::queue<int, int, int[]>(world);
	for (auto i = 0u; i < indices_samples.size(); i++) {
		for (int t = 0; t < p; t++) {
			sample_queue(t).send(s, i, H(indices_samples[i]).nets());
		}
	}
	world.sync();
	
	int total_samples = p * opts.sample_size;
	
	/* Compute the inner products of the samples and the local vertices */
	auto ip = std::vector<std::vector<double>>(H.size(), std::vector<double>(total_samples, 0.0));
	auto degree_samples = std::vector<int>(total_samples);

	for (auto& [t, number_sample, sample_nets] : std::move(sample_queue)) {
		degree_samples[t * opts.sample_size + number_sample] = sample_nets.size();
		for (auto n_id : sample_nets) {
			for (auto u_id : H.net(n_id).vertices()) {
				ip[H.local_id(u_id)][t * opts.sample_size + number_sample] += 1.0/((double)H.net(n_id).global_size() - 1.0);
			}
		}
	}
	
	/* We set the ip of all local samples with itself to 0, so they will not match themselves */
	for (auto i = 0u; i < indices_samples.size(); i++) {
		ip[indices_samples[i]][s * opts.sample_size + i] = 0;
	}
	
	/* Find best sample for vertex v and add it to the list of that sample */
	auto requested_matches = std::vector<std::vector<std::pair<int, double>>>(total_samples);
	for (auto& v : H.vertices()) {
		double max_ip = 0.0;
		int best_match = -1;
		auto local_id = H.local_id(v.id());
		for (auto u = 0; u < total_samples; u++) {
			double ip_vu = ip[local_id][u];
			if (ip_vu > 0) {
				ip_vu *= 1.0/(double)std::min((int)v.degree(), degree_samples[u]);
				if (ip_vu > max_ip) {
					max_ip = ip_vu;
					best_match = u;
				}
			}
		}
		if (best_match != -1) {
			requested_matches[best_match].push_back(std::make_pair(local_id, max_ip));
		}
	}
	
	for (auto& match_list : requested_matches) {
		std::sort(match_list.begin(), match_list.end(), [](const auto& match1, const auto& match2) -> bool { return match1.second > match2.second; });
	}
	
	/* Queue for the vertex requests with the vertex to match with, the id of the vertex that wants to match and their ip */
	auto request_queue = bulk::queue<int, int, double>(world);
	for (int sample = 0; sample < total_samples; sample++) {
		int t = sample/opts.sample_size;
		int number_to_send = std::min((int)requested_matches[sample].size(), opts.coarsening_max_clustersize);
		for (int i = 0; i < number_to_send; i++) {
			request_queue(t).send(sample - t * opts.sample_size, H(requested_matches[sample][i].first).id(), requested_matches[sample][i].second);
		}
	}
	
	world.sync();
	
	for (auto& [sample, proposer, scip] : request_queue) {
		world.log("sample global id: %d, proposer: %d, scip: %lf", H(indices_samples[sample]).id(), proposer, scip);
	}
	
	world.sync();

	return H;
}

} // namespace pmondriaan
