#include <algorithm>
#include <array>
#include <iostream>
#include <unordered_set>
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

    // queue to send the information about the accepted samples
    auto info_queue = bulk::queue<int, long, int[], long[]>(world);
    auto matched = std::vector<bool>(H.size(), false);
    pmondriaan::send_information_matches(H, accepted_matches, info_queue,
                                         matched, opts.sample_size);

    auto HC = pmondriaan::contract_hypergraph(world, H, indices_samples, info_queue, matched);
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

/**
 * First merges the nets and weight of all vertices matched to a sample and
 * then sends this information to the owner of the sample.
 */
void send_information_matches(pmondriaan::hypergraph& H,
                              bulk::queue<int, int>& accepted_matches,
                              bulk::queue<int, long, int[], long[]>& info_queue,
                              std::vector<bool>& matched,
                              int sample_size) {

    std::sort(accepted_matches.begin(), accepted_matches.end());
    int prev_sample = -1;
    long total_weight_sample = 0;
    auto total_nets_sample = std::unordered_set<int>();
    for (auto& [sample, proposer] : accepted_matches) {
        if ((sample != prev_sample) && (prev_sample >= 0)) {
            // We send all information about the matches of the previous sample to the correct processor
            int t = prev_sample / sample_size;
            auto nets_vector = std::vector<int>();
            auto cost_nets = std::vector<long>();
            for (auto n : total_nets_sample) {
                nets_vector.push_back(n);
                cost_nets.push_back(H.net(n).cost());
            }

            info_queue(t).send(sample - t * sample_size, total_weight_sample,
                               nets_vector, cost_nets);
            total_weight_sample = 0;
            total_nets_sample.clear();
            prev_sample = sample;
        }
        matched[H.local_id(proposer)] = true;
        auto v = H(H.local_id(proposer));
        total_weight_sample += v.weight();
        total_nets_sample.insert(v.nets().begin(), v.nets().end());
    }
    info_queue.world().sync();
}

pmondriaan::hypergraph contract_hypergraph(bulk::world& world,
                                           pmondriaan::hypergraph& H,
                                           const std::vector<int> samples,
                                           bulk::queue<int, long, int[], long[]>& matches,
                                           std::vector<bool>& matched) {

    // we new nets to which we will later add the vertices
    auto new_nets = std::vector<pmondriaan::net>();
    std::unordered_map<int, int> net_global_to_local;
    size_t number_nets = 0;

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
                auto insert_result = net_global_to_local.insert({n, number_nets});
                if (insert_result.second) {
                    new_nets.push_back(
                    pmondriaan::net(n, std::vector<int>(), H.net(n).cost()));
                    number_nets++;
                }
                new_nets[net_global_to_local[n]].add_vertex(v.id());
            }
        }
    }

    // we create the new weight and adjacency list for all samples
    auto sample_total_weight = std::vector<long>(samples.size(), 0);
    auto sample_net_lists = std::vector<std::unordered_set<int>>(samples.size());
    for (const auto& [sample, weight, nets, cost_nets] : matches) {
        sample_total_weight[sample] += weight;
        sample_net_lists[sample].insert(nets.begin(), nets.end());
        for (auto i = 0u; i < nets.size(); i++) {
            auto insert_result = net_global_to_local.insert({nets[i], number_nets});
            if (insert_result.second) {
                new_nets.push_back(
                pmondriaan::net(nets[i], std::vector<int>(), cost_nets[i]));
                number_nets++;
            }
        }
    }

    int new_size = (int)new_vertices.size();
    auto new_global_size = bulk::sum(matches.world(), new_size);
    auto HC = pmondriaan::hypergraph(new_global_size, H.global_number_nets(),
                                     new_vertices, new_nets);


    for (auto index = 0u; index < samples.size(); index++) {
        int i = (int)index;
        auto& sample_vertex = HC(HC.local_id(H(samples[i]).id()));
        auto& nets = sample_vertex.nets();
        std::unordered_set<int> nets_set(nets.begin(), nets.end());

        sample_vertex.add_weight(sample_total_weight[i]);

        for (auto net : sample_net_lists[i]) {
            auto insert_result = nets_set.insert(net);
            if (insert_result.second) {
                nets.push_back(net);
                HC.net(net).add_vertex(sample_vertex.id());
            }
        }
    }

    pmondriaan::global_net_sizes(world, HC);
    remove_free_nets(world, HC);

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
    world.sync();
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

    int sample_count = 0;
    for (auto& new_vertex : HC.vertices()) {
        auto local_id = H.local_id(new_vertex.id());

        C.add_sample(new_vertex.id());
        sample_count++;

        auto& nets = new_vertex.nets();
        for (auto n : nets) {
            HC.net(n).add_vertex(new_vertex.id());
        }

        if (matches[local_id].size() > 0) {
            std::unordered_set<int> nets_set(nets.begin(), nets.end());

            for (auto match : matches[local_id]) {
                C.add_match(sample_count - 1, match, world.rank());
                new_vertex.add_weight(H(H.local_id(match)).weight());
                for (auto net : H(H.local_id(match)).nets()) {
                    auto insert_result = nets_set.insert(net);
                    if (insert_result.second) {
                        nets.push_back(net);
                        HC.net(net).add_vertex(new_vertex.id());
                    }
                }
            }
        }
    }

    remove_free_nets(HC);

    for (auto& net : HC.nets()) {
        net.set_global_size(net.size());
    }

    return HC;
} // namespace pmondriaan

} // namespace pmondriaan
