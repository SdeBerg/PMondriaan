#include <vector>
#include <array>
#include <algorithm>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "multilevel_bisect/coarsen.hpp"
#include "bisect.hpp"
#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/sample.hpp"
#include "hypergraph/contraction.hpp"

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
	
	world.log("samples: %d", indices_samples.size());
	
	/* We now send the samples and the processor id to all processors */
	auto sample_queue = bulk::queue<int, long, int[]>(world);
	for (auto i = 0u; i < indices_samples.size(); i++) {
		for (int t = 0; t < p; t++) {
			sample_queue(t).send(s, (long)i, H(indices_samples[i]).nets());
		}
	}

	world.sync();
	
	auto accepted_matches = bulk::queue<int, int>(world);
	/* After his funtion, accepted matches contains the matches that have been accepted */
	request_matches(H, sample_queue, accepted_matches, indices_samples, opts);

	/* We use the sample_queue again to send the information about the accepted samples */
	auto matched = std::vector<bool>(H.size(), false);
	for (auto& [sample, proposer] : accepted_matches) {
		matched[H.local_id(proposer)] = true;
		auto v = H(H.local_id(proposer));
		int t = sample / opts.sample_size;
		sample_queue(t).send(sample - t * opts.sample_size, v.weight(), v.nets());
	}
	
	world.sync();
	
	/* We compute the number of new vertices (in the matched vector samples are also labeled as unmatched */
	int new_vertices = 0;
	for (auto m : matched) {
		if (!m) {
			new_vertices++;
		}
	}
	//auto contraction = pmondriaan::contraction(new_vertices);
	auto HC = pmondriaan::contract_hypergraph(world, H, indices_samples, sample_queue, matched);
	return HC;
}


/* Sends match request to the owners of the best matches found using the improduct computation.
 * Returns the local matches. 
 */
void request_matches(pmondriaan::hypergraph& H, auto& sample_queue, auto& accepted_matches, const std::vector<int>& indices_samples, pmondriaan::options& opts) {
	
	auto& world = sample_queue.world();
	int s = world.rank();
	int p = world.active_processors();
	size_t number_local_samples = indices_samples.size();
	
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
	
	/* We set the ip of all local samples with all samples to 0, so they will not match eachother */
	for (auto i = 0u; i < number_local_samples; i++) {
		for (auto j = 0; j < total_samples; j++) {
			ip[indices_samples[i]][j] = 0.0;
		}
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
			requested_matches[best_match].push_back(std::make_pair(v.id(), max_ip));
		}
	}
	
	for (auto& match_list : requested_matches) {
		std::sort(match_list.begin(), match_list.end(), [](const auto& match1, const auto& match2) -> bool { return match1.second > match2.second; });
	}
	
	/* Queue for the vertex requests with the sender, the vertex to match with, the id of the vertex that wants to match and their ip */
	auto request_queue = bulk::queue<int, int, int, double>(world);
	for (int sample = 0; sample < total_samples; sample++) {
		int t = sample/opts.sample_size;
		int number_to_send = std::min((int)requested_matches[sample].size(), opts.coarsening_max_clustersize);
		for (int i = 0; i < number_to_send; i++) {
			request_queue(t).send(s, sample - t * opts.sample_size, requested_matches[sample][i].first, requested_matches[sample][i].second);
		}
	}
	
	world.sync();
	
	auto matches = std::vector<std::vector<std::tuple<int, int, double>>>(number_local_samples);
	for (const auto& [sender, sample, proposer, scip] : request_queue) {
		matches[sample].push_back(std::make_tuple(sender, proposer, scip));
	}
	
	for (auto i = 0u; i < indices_samples.size(); i++) {
		auto& match_list = matches[i];
		std::sort(match_list.begin(), match_list.end(), [](const auto& match1, const auto& match2) -> bool { return std::get<2>(match1) > std::get<2>(match2); });
		
		int number_to_send = std::min((int)match_list.size(), opts.coarsening_max_clustersize);
		for (int j = 0; j < number_to_send; j++) {
			auto& match = match_list[j];
			accepted_matches(std::get<0>(match)).send(i + s * opts.sample_size, std::get<1>(match));
		}
	}
	
	world.sync();
}

pmondriaan::hypergraph contract_hypergraph(bulk::world& world, pmondriaan::hypergraph& H, const std::vector<int> samples, auto& matches, std::vector<bool>& matched) {

	/* We new nets to which we will later add the vertices */
	auto new_nets = std::vector<pmondriaan::net>();
	for (auto& net : H.nets()) {
		new_nets.push_back(pmondriaan::net(net.id(), std::vector<int>(), net.cost()));
	}
	/* For each unmatched vertex we create a new one (samples are also "unmatched") */
	auto new_vertices = std::vector<pmondriaan::vertex>();
	for (auto index = 0u; index < H.size(); index++) {
		if (!matched[index]) {
			auto& v = H(index);
			
			auto new_v_nets = std::vector<int>();
			new_v_nets.insert(new_v_nets.begin(), v.nets().begin(), v.nets().end());
			auto new_v = pmondriaan::vertex(v.id(), new_v_nets, v.weight());
			new_vertices.push_back(new_v);
			
			for (auto n : new_v.nets()) {
				new_nets[n].add_vertex(v.id());
			}
		}
	}
	
	bulk::var<int> new_size(matches.world());
	new_size = (int)new_vertices.size();
	auto new_global_size = bulk::foldl(new_size, [](int& lhs, int rhs) { lhs += rhs; });
	
	auto HC = pmondriaan::hypergraph(new_global_size, new_vertices, new_nets);
	
	/* We create the new weight and adjacency list for all samples */
	auto sample_total_weight= std::vector<long>(samples.size(), 0);
	auto sample_net_lists = std::vector<std::vector<int>>(samples.size());
	for (auto& [sample, weight, nets] : std::move(matches)) {
		sample_total_weight[sample] += weight;
		sample_net_lists[sample].insert(sample_net_lists[sample].end(), nets.begin(), nets.end());
	}
	
	/* This vector keeps track of the nets already include for the current vertex
	 if included_in_net[n] is equal to the sample number, the net is already included */
	auto included_in_net = std::vector<int>(H.nets().size(), -1);
	for (auto index = 0u; index < samples.size(); index++) {
		int i = (int)index;
		auto& sample_vertex = HC(HC.local_id(H(samples[i]).id()));
		auto& nets = sample_vertex.nets();
		
		sample_vertex.add_weight(sample_total_weight[i]);
		
		for (auto n : nets) {
			included_in_net[n] = i;
		}
		for (auto net : sample_net_lists[i]) {
			if (included_in_net[net] != i) {
				nets.push_back(net);
				HC.net(net).add_vertex(sample_vertex.id());
				included_in_net[net] = i;
			}
		}
	}
	
	pmondriaan::global_net_sizes(world, HC);
	
	return HC;
}








} // namespace pmondriaan
