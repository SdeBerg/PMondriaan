#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "bisect.hpp"
#include "hypergraph/contraction.hpp"
#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/coarsen.hpp"
#include "multilevel_bisect/sample.hpp"

namespace pmondriaan {

/**
 * Coarsens the hypergraph H and returns a new hypergraph HC.
 */
pmondriaan::hypergraph coarsen_hypergraph_par(bulk::world& world,
                                              pmondriaan::hypergraph& H,
                                              pmondriaan::contraction& C,
                                              pmondriaan::options& opts,
                                              std::mt19937& rng) {

    int s = world.rank();
    int p = world.active_processors();

    /* We first select ns samples */
    auto indices_samples = std::vector<int>();
    if (opts.sampling_mode == pmondriaan::sampling::random) {
        indices_samples = sample_random(H, opts.sample_size, rng);
    } else if (opts.sampling_mode == pmondriaan::sampling::label_propagation) {
        indices_samples = sample_lp(H, opts, rng);
    }

    world.log("samples: %d", indices_samples.size());

    // we now send the samples and the processor id to all processors
    auto sample_queue = bulk::queue<int, long, int[]>(world);
    for (auto i = 0u; i < indices_samples.size(); i++) {
        for (int t = 0; t < p; t++) {
            sample_queue(t).send(s, (long)i, H(indices_samples[i]).nets());
        }
    }

    world.sync();

    C.add_samples(H, indices_samples);
    auto accepted_matches = bulk::queue<int, int>(world);
    // after his funtion, accepted matches contains the matches that have been accepted
    request_matches(H, C, sample_queue, accepted_matches, indices_samples, opts);

    // we use the sample_queue again to send the information about the accepted samples
    auto matched = std::vector<bool>(H.size(), false);
    for (auto& [sample, proposer] : accepted_matches) {
        matched[H.local_id(proposer)] = true;
        auto v = H(H.local_id(proposer));
        int t = sample / opts.sample_size;
        sample_queue(t).send(sample - t * opts.sample_size, v.weight(), v.nets());
    }

    world.sync();

    auto HC =
    pmondriaan::contract_hypergraph(world, H, indices_samples, sample_queue, matched);
    return HC;
}


/**
 * Sends match request to the owners of the best matches found using the
 * improduct computation. Returns the local matches.
 */
void request_matches(pmondriaan::hypergraph& H,
                     pmondriaan::contraction& C,
                     bulk::queue<int, long, int[]>& sample_queue,
                     bulk::queue<int, int>& accepted_matches,
                     const std::vector<int>& indices_samples,
                     pmondriaan::options& opts) {

    auto& world = sample_queue.world();
    int s = world.rank();
    int p = world.active_processors();
    size_t number_local_samples = indices_samples.size();

    int total_samples = p * opts.sample_size;

    // compute the inner products of the samples and the local vertices
    auto ip =
    std::vector<std::vector<double>>(H.size(), std::vector<double>(total_samples, 0.0));
    auto degree_samples = std::vector<int>(total_samples);

    for (const auto& [t, number_sample, sample_nets] : sample_queue) {
        degree_samples[t * opts.sample_size + number_sample] = sample_nets.size();
        for (auto n_id : sample_nets) {
            double scaled_cost = H.net(n_id).scaled_cost();
            for (auto u_id : H.net(n_id).vertices()) {
                ip[H.local_id(u_id)][t * opts.sample_size + number_sample] += scaled_cost;
            }
        }
    }

    // we set the ip of all local samples with all samples to 0, so they will not match eachother
    for (auto i = 0u; i < number_local_samples; i++) {
        for (auto j = 0; j < total_samples; j++) {
            ip[indices_samples[i]][j] = 0.0;
        }
    }

    // find best sample for vertex v and add it to the list of that sample
    auto requested_matches = std::vector<std::vector<std::pair<int, double>>>(total_samples);
    for (auto& v : H.vertices()) {
        double max_ip = 0.0;
        int best_match = -1;
        auto local_id = H.local_id(v.id());
        for (auto u = 0; u < total_samples; u++) {
            double ip_vu = ip[local_id][u];
            if (ip_vu > 0) {
                ip_vu *= (1.0 / (double)std::min((int)v.degree(), degree_samples[u]));
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
        std::sort(match_list.begin(), match_list.end(),
                  [](const auto& match1, const auto& match2) -> bool {
                      return match1.second > match2.second;
                  });
    }

    // queue for the vertex requests with the sender, the vertex to match with, the id of the vertex that wants to match and their ip
    auto request_queue = bulk::queue<int, int, int, double>(world);
    for (int sample = 0; sample < total_samples; sample++) {
        int t = sample / opts.sample_size;
        int number_to_send = std::min((int)requested_matches[sample].size(),
                                      opts.coarsening_max_clustersize);
        for (int i = 0; i < number_to_send; i++) {
            request_queue(t).send(s, sample - t * opts.sample_size,
                                  requested_matches[sample][i].first,
                                  requested_matches[sample][i].second);
        }
    }

    world.sync();

    auto matches =
    std::vector<std::vector<std::tuple<int, int, double>>>(number_local_samples);
    for (const auto& [sender, sample, proposer, scip] : request_queue) {
        matches[sample].push_back(std::make_tuple(sender, proposer, scip));
    }

    for (auto i = 0u; i < number_local_samples; i++) {
        auto& match_list = matches[i];
        std::sort(match_list.begin(), match_list.end(),
                  [](const auto& match1, const auto& match2) -> bool {
                      return std::get<2>(match1) > std::get<2>(match2);
                  });

        int number_to_send =
        std::min((int)match_list.size(), opts.coarsening_max_clustersize);
        for (int j = 0; j < number_to_send; j++) {
            auto& match = match_list[j];
            C.add_match(i, std::get<1>(match), std::get<0>(match));
            accepted_matches(std::get<0>(match)).send(i + s * opts.sample_size, std::get<1>(match));
        }
    }

    world.sync();
}

pmondriaan::hypergraph contract_hypergraph(bulk::world& world,
                                           pmondriaan::hypergraph& H,
                                           const std::vector<int> samples,
                                           bulk::queue<int, long, int[]>& matches,
                                           std::vector<bool>& matched) {

    // we new nets to which we will later add the vertices
    auto new_nets = std::vector<pmondriaan::net>();
    for (auto& net : H.nets()) {
        new_nets.push_back(pmondriaan::net(net.id(), std::vector<int>(), net.cost()));
    }
    // for each unmatched vertex we create a new one (samples are also "unmatched")
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

    // we create the new weight and adjacency list for all samples
    auto sample_total_weight = std::vector<long>(samples.size(), 0);
    auto sample_net_lists = std::vector<std::vector<int>>(samples.size());
    for (auto& [sample, weight, nets] : std::move(matches)) {
        sample_total_weight[sample] += weight;
        sample_net_lists[sample].insert(sample_net_lists[sample].end(),
                                        nets.begin(), nets.end());
    }

    int new_size = (int)new_vertices.size();
    auto new_global_size = bulk::sum(matches.world(), new_size);
    auto HC = pmondriaan::hypergraph(new_global_size, H.global_number_nets(),
                                     new_vertices, new_nets);

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


/**
 * Coarsens the hypergraph H and returns a hypergraph HC sequentially.
 */
pmondriaan::hypergraph coarsen_hypergraph_seq(bulk::world& world,
                                              pmondriaan::hypergraph& H,
                                              pmondriaan::contraction& C,
                                              pmondriaan::options& opts,
                                              std::mt19937& rng) {

    auto matches = std::vector<std::vector<int>>(H.size(), std::vector<int>());
    auto matched = std::vector<bool>(H.size(), false);
    // contains the vertices of the contracted hypergraph
    auto new_v = std::vector<pmondriaan::vertex>();

    auto ip = std::vector<double>(H.size(), 0.0);
    // we visit the vertices in a random order
    std::vector<int> indices(H.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    for (auto i : indices) {
        auto& v = H(i);
        if (matches[i].empty()) {
            auto visited = std::vector<int>();
            for (auto n_id : v.nets()) {
                double scaled_cost = H.net(n_id).scaled_cost();
                for (auto u_id : H.net(n_id).vertices()) {
                    auto u_local = H.local_id(u_id);
                    if ((!matched[u_local]) && (u_local != i)) {
                        if (ip[u_local] == 0.0) {
                            visited.push_back(u_local);
                        }
                        ip[u_local] += scaled_cost;
                    }
                }
            }

            double max_ip = 0.0;
            int best_match = -1;
            for (auto u : visited) {
                ip[u] *= (1.0 / (double)std::min(v.degree(), H(u).degree()));
                if ((ip[u] > max_ip) &&
                    ((int)matches[u].size() < opts.coarsening_max_clustersize)) {
                    max_ip = ip[u];
                    best_match = u;
                }
            }
            if (best_match != -1) {
                matches[best_match].push_back(v.id());
                matched[i] = true;
            } else {
                add_v_to_list(new_v, v);
            }

            for (auto u : visited) {
                ip[u] = 0.0;
            }
        } else {
            add_v_to_list(new_v, v);
        }
    }
    auto result = pmondriaan::contract_hypergraph(world, H, C, matches, new_v);
    return result;
}

void add_v_to_list(std::vector<pmondriaan::vertex>& v_list, pmondriaan::vertex& v) {
    auto new_v_nets = std::vector<int>();
    new_v_nets.insert(new_v_nets.begin(), v.nets().begin(), v.nets().end());
    v_list.push_back(pmondriaan::vertex(v.id(), new_v_nets, v.weight()));
}

pmondriaan::hypergraph contract_hypergraph(bulk::world& world,
                                           pmondriaan::hypergraph& H,
                                           pmondriaan::contraction& C,
                                           std::vector<std::vector<int>>& matches,
                                           std::vector<pmondriaan::vertex>& new_vertices) {

    // we new nets to which we will later add the vertices
    auto new_nets = std::vector<pmondriaan::net>();
    for (auto& net : H.nets()) {
        new_nets.push_back(pmondriaan::net(net.id(), std::vector<int>(), net.cost()));
    }

    auto HC = pmondriaan::hypergraph(new_vertices.size(), H.global_number_nets(),
                                     new_vertices, new_nets);

    // This vector keeps track of the nets already include for the current vertex
    // if included_in_net[n] is equal to the sample number, the net is already included
    auto included_in_net = std::vector<int>(H.nets().size(), -1);
    int sample_count = 0;
    for (auto index = 0u; index < H.size(); index++) {
        int i = (int)index;
        if (matches[i].size() > 0) {
            C.add_sample(H(i).id());
            sample_count++;

            auto& nets = HC(HC.local_id(H(i).id())).nets();
            for (auto n : nets) {
                included_in_net[n] = i;
            }

            for (auto match : matches[i]) {
                C.add_match(sample_count - 1, match, world.rank());
                HC(HC.local_id(H(i).id())).add_weight(H(H.local_id(match)).weight());
                for (auto net : H(H.local_id(match)).nets()) {
                    if (included_in_net[net] != i) {
                        nets.push_back(net);
                        HC.net(net).add_vertex(H(i).id());
                        included_in_net[net] = i;
                    }
                }
            }
        }
    }

    for (auto& net : HC.nets()) {
        net.set_global_size(net.size());
    }

    return HC;
}

} // namespace pmondriaan
