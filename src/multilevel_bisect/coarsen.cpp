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
    auto s = world.rank();
    auto p = world.active_processors();

    /* We first select ns samples */
    auto indices_samples = std::vector<long>();
    if (opts.sampling_mode == pmondriaan::sampling::random) {
        indices_samples = sample_random(H, opts.sample_size, rng);
    } else if (opts.sampling_mode == pmondriaan::sampling::label_propagation) {
        indices_samples = sample_lp(H, opts, rng);
    }

    // we now send the samples and the processor id to all processors
    auto sample_queue = bulk::queue<long, long, long[]>(world);
    for (auto i = 0u; i < indices_samples.size(); i++) {
        for (long t = 0; t < p; t++) {
            sample_queue(t).send(s, (long)i, H(indices_samples[i]).nets());
        }
    }

    world.sync();

    C.add_samples(H, indices_samples);
    auto accepted_matches = bulk::queue<long, long>(world);
    // after his funtion, accepted matches contains the matches that have been accepted

    request_matches(H, C, sample_queue, accepted_matches, indices_samples, opts);

    // queue to send the information about the accepted samples
    auto info_queue = bulk::queue<long, long, long[], long[]>(world);
    auto matched = std::vector<bool>(H.size(), false);

    pmondriaan::send_information_matches(world, H, accepted_matches, info_queue,
                                         matched, opts.sample_size);

    auto HC =
    pmondriaan::contract_hypergraph(world, H, C, indices_samples, info_queue, matched);

    return HC;
}


/**
 * Sends match request to the owners of the best matches found using the
 * improduct computation. Returns the local matches.
 */
void request_matches(pmondriaan::hypergraph& H,
                     pmondriaan::contraction& C,
                     bulk::queue<long, long, long[]>& sample_queue,
                     bulk::queue<long, long>& accepted_matches,
                     const std::vector<long>& indices_samples,
                     pmondriaan::options opts) {

    auto& world = sample_queue.world();
    auto s = world.rank();
    auto p = world.active_processors();
    size_t number_local_samples = indices_samples.size();

    long total_samples = p * opts.sample_size;

    // compute the inner products of the samples and the local vertices
    auto best_ip =
    std::vector<std::pair<double, long>>(H.size(), std::make_pair(0.0, -1));
    auto current_ip = std::vector<double>(H.size(), 0.0);
    std::unordered_set<long> changed_indices;

    for (const auto& [t, number_sample, sample_nets] : sample_queue) {
        for (auto n_id : sample_nets) {
            if (H.is_local_net(n_id)) {
                double scaled_cost = H.net(n_id).scaled_cost();
                for (auto u_id : H.net(n_id).vertices()) {
                    assert(H.local_id(u_id) >= 0 &&
                           (size_t)H.local_id(u_id) < current_ip.size());
                    current_ip[H.local_id(u_id)] += scaled_cost;
                    changed_indices.insert(H.local_id(u_id));
                }
            }
        }
        for (auto index : changed_indices) {
            assert(index >= 0 && (size_t)index < current_ip.size() &&
                   (size_t)index < best_ip.size());
            size_t min_degree = std::min(H(index).degree(), sample_nets.size());
            assert(min_degree > 0);
            current_ip[index] *= 1.0 / (double)min_degree;
            if (current_ip[index] > best_ip[index].first) {
                best_ip[index] =
                std::make_pair(current_ip[index], t * opts.sample_size + number_sample);
            }
            current_ip[index] = 0.0;
        }
        changed_indices.clear();
    }

    // we set the ip of all local samples with all samples to 0, so they will not match eachother
    for (auto local_sample : indices_samples) {
        assert(local_sample >= 0 && (size_t)local_sample < best_ip.size());
        best_ip[local_sample] = std::make_pair(0.0, -1);
    }

    // find best sample for vertex v and add it to the list of that sample
    auto requested_matches =
    std::vector<std::vector<std::pair<long, double>>>(total_samples);
    for (auto& v : H.vertices()) {
        auto local_id = H.local_id(v.id());
        if (best_ip[local_id].second != -1) {
            assert(best_ip[local_id].second >= 0 &&
                   (size_t)best_ip[local_id].second < requested_matches.size());

            assert(local_id >= 0 && (size_t)local_id < best_ip.size());

            requested_matches[best_ip[local_id].second].push_back(
            std::make_pair(v.id(), best_ip[local_id].first));
        }
    }

    for (auto& match_list : requested_matches) {
        if (match_list.size() > opts.coarsening_max_clustersize) {
            assert(match_list.begin() + opts.coarsening_max_clustersize <
                   match_list.end());
            std::nth_element(match_list.begin(), match_list.begin() + opts.coarsening_max_clustersize,
                             match_list.end(),
                             [](const auto& match1, const auto& match2) -> bool {
                                 return match1.second > match2.second;
                             });
        }
    }

    // queue for the vertex requests with the sender, the vertex to match with, the id of the vertex that wants to match and their ip
    auto request_queue = bulk::queue<long, long, long, double>(world);
    for (long sample = 0; sample < total_samples; sample++) {
        assert(opts.sample_size > 0);
        long t = sample / opts.sample_size;
        assert(t >= 0 && t < p);
        assert(sample >= 0 && (size_t)sample < requested_matches.size());
        long number_to_send =
        std::min(requested_matches[sample].size(), opts.coarsening_max_clustersize);
        for (long i = 0; i < number_to_send; i++) {
            request_queue(t).send(s, sample - t * opts.sample_size,
                                  requested_matches[sample][i].first,
                                  requested_matches[sample][i].second);
        }
    }

    world.sync();

    auto matches =
    std::vector<std::vector<std::tuple<long, long, double>>>(number_local_samples);
    for (const auto& [sender, sample, proposer, scip] : request_queue) {
        assert((size_t)sample < number_local_samples);
        matches[sample].push_back(std::make_tuple(sender, proposer, scip));
    }

    for (auto i = 0u; i < number_local_samples; i++) {
        auto& match_list = matches[i];
        if (match_list.size() > opts.coarsening_max_clustersize) {
            assert(match_list.begin() + opts.coarsening_max_clustersize <
                   match_list.end());
            std::nth_element(match_list.begin(), match_list.begin() + opts.coarsening_max_clustersize,
                             match_list.end(),
                             [](const auto& match1, const auto& match2) -> bool {
                                 return std::get<2>(match1) > std::get<2>(match2);
                             });
        }

        long number_to_send = std::min(match_list.size(), opts.coarsening_max_clustersize);
        for (long j = 0; j < number_to_send; j++) {
            auto& match = match_list[j];
            C.add_match(i, std::get<1>(match), std::get<0>(match));
            auto t = std::get<0>(match);
            assert(t >= 0 && t < p);
            accepted_matches(std::get<0>(match)).send(i + s * opts.sample_size, std::get<1>(match));
        }
    }
    world.sync();
}

/**
 * First merges the nets and weight of all vertices matched to a sample and
 * then sends this information to the owner of the sample.
 */
void send_information_matches(bulk::world& world,
                              pmondriaan::hypergraph& H,
                              bulk::queue<long, long>& accepted_matches,
                              bulk::queue<long, long, long[], long[]>& info_queue,
                              std::vector<bool>& matched,
                              long sample_size) {
    std::sort(accepted_matches.begin(), accepted_matches.end());
    long prev_sample = -1;
    long total_weight_sample = 0;
    auto total_nets_sample = std::unordered_set<long>();
    for (auto& [sample, proposer] : accepted_matches) {
        if ((sample != prev_sample) && (prev_sample >= 0)) {
            assert(sample_size > 0);
            // We send all information about the matches of the previous sample to the correct processor
            long t = prev_sample / sample_size;
            auto nets_vector = std::vector<long>();
            auto cost_nets = std::vector<long>();
            for (auto n : total_nets_sample) {
                nets_vector.push_back(n);
                cost_nets.push_back(H.net(n).cost());
            }
            assert(t >= 0 && t < world.active_processors());
            info_queue(t).send(prev_sample - t * sample_size,
                               total_weight_sample, nets_vector, cost_nets);
            total_weight_sample = 0;
            total_nets_sample.clear();
        }
        assert(H.local_id(proposer) >= 0 &&
               (size_t)H.local_id(proposer) < matched.size());
        matched[H.local_id(proposer)] = true;
        auto v = H(H.local_id(proposer));
        total_weight_sample += v.weight();
        total_nets_sample.insert(v.nets().begin(), v.nets().end());
        prev_sample = sample;
    }

    // We send all information about the last sample
    assert(sample_size > 0);
    long t = prev_sample / sample_size;
    auto nets_vector = std::vector<long>();
    auto cost_nets = std::vector<long>();
    for (auto n : total_nets_sample) {
        nets_vector.push_back(n);
        cost_nets.push_back(H.net(n).cost());
    }
    assert(t >= 0 && t < world.active_processors());
    info_queue(t).send(prev_sample - t * sample_size, total_weight_sample,
                       nets_vector, cost_nets);

    world.sync();
}

pmondriaan::hypergraph contract_hypergraph(bulk::world& world,
                                           pmondriaan::hypergraph& H,
                                           pmondriaan::contraction& C,
                                           const std::vector<long> samples,
                                           bulk::queue<long, long, long[], long[]>& matches,
                                           std::vector<bool>& matched) {

    // new nets to which we will later add the vertices
    auto new_nets = std::vector<pmondriaan::net>();
    std::unordered_map<long, long> net_global_to_local;
    size_t number_nets = 0;

    // for each unmatched vertex we create a new one (samples are also "unmatched")
    auto new_vertices = std::vector<pmondriaan::vertex>();
    for (auto index = 0u; index < H.size(); index++) {
        if (!matched[index]) {
            auto& v = H(index);

            auto new_v_nets = std::vector<long>();
            new_v_nets.insert(new_v_nets.begin(), v.nets().begin(), v.nets().end());
            auto new_v = pmondriaan::vertex(v.id(), new_v_nets, v.weight());
            new_vertices.push_back(new_v);

            for (auto n : new_v.nets()) {
                auto insert_result = net_global_to_local.insert({n, number_nets});
                if (insert_result.second) {
                    new_nets.push_back(
                    pmondriaan::net(n, std::vector<long>(), H.net(n).cost()));
                    number_nets++;
                }
                new_nets[net_global_to_local[n]].add_vertex(v.id());
            }
        }
    }

    // we create the new weight and adjacency list for all samples
    auto sample_total_weight = std::vector<long>(samples.size(), 0);
    auto sample_net_lists = std::vector<std::unordered_set<long>>(samples.size());
    for (const auto& [sample, weight, nets, cost_nets] : matches) {
        sample_total_weight[sample] += weight;
        sample_net_lists[sample].insert(nets.begin(), nets.end());
        for (auto i = 0u; i < nets.size(); i++) {
            auto insert_result = net_global_to_local.insert({nets[i], number_nets});
            if (insert_result.second) {
                new_nets.push_back(
                pmondriaan::net(nets[i], std::vector<long>(), cost_nets[i]));
                number_nets++;
            }
        }
    }

    auto new_size = new_vertices.size();
    auto new_global_size = bulk::sum(matches.world(), new_size);
    auto HC = pmondriaan::hypergraph(new_global_size, H.global_number_nets(),
                                     new_vertices, new_nets);


    for (auto index = 0u; index < samples.size(); index++) {
        long i = (long)index;
        auto& sample_vertex = HC(HC.local_id(H(samples[i]).id()));
        auto& nets = sample_vertex.nets();
        std::unordered_set<long> nets_set(nets.begin(), nets.end());

        sample_vertex.add_weight(sample_total_weight[i]);

        for (auto net : sample_net_lists[i]) {
            auto insert_result = nets_set.insert(net);
            if (insert_result.second) {
                nets.push_back(net);
                HC.net(net).add_vertex(sample_vertex.id());
            }
        }
    }

    remove_free_nets(world, HC, 1);
    C.merge_free_vertices(world, HC);

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
    auto matches = std::vector<std::vector<long>>(H.size(), std::vector<long>());
    auto matched = std::vector<bool>(H.size(), false);
    // contains the vertices of the contracted hypergraph
    auto new_v = std::vector<pmondriaan::vertex>();

    auto ip = std::vector<double>(H.size(), 0.0);
    // we visit the vertices in a random order
    std::vector<long> indices(H.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    for (auto i : indices) {
        auto& v = H(i);
        if (matches[i].empty()) {
            auto visited = std::vector<long>();
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
            long best_match = -1;
            for (auto u : visited) {
                ip[u] *= (1.0 / (double)std::min(v.degree(), H(u).degree()));
                if ((ip[u] > max_ip) && (matches[u].size() < opts.coarsening_max_clustersize)) {
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
    auto new_v_nets = std::vector<long>();
    new_v_nets.insert(new_v_nets.begin(), v.nets().begin(), v.nets().end());
    v_list.push_back(pmondriaan::vertex(v.id(), new_v_nets, v.weight()));
}

pmondriaan::hypergraph contract_hypergraph(bulk::world& world,
                                           pmondriaan::hypergraph& H,
                                           pmondriaan::contraction& C,
                                           std::vector<std::vector<long>>& matches,
                                           std::vector<pmondriaan::vertex>& new_vertices) {
    // we new nets to which we will later add the vertices
    auto new_nets = std::vector<pmondriaan::net>();
    for (auto& net : H.nets()) {
        new_nets.push_back(pmondriaan::net(net.id(), std::vector<long>(), net.cost()));
    }
    auto HC = pmondriaan::hypergraph(new_vertices.size(), H.global_number_nets(),
                                     new_vertices, new_nets);

    long sample_count = 0;
    for (auto& new_vertex : HC.vertices()) {
        auto local_id = H.local_id(new_vertex.id());

        C.add_sample(new_vertex.id());
        sample_count++;

        auto& nets = new_vertex.nets();
        for (auto n : nets) {
            HC.net(n).add_vertex(new_vertex.id());
        }

        if (matches[local_id].size() > 0) {
            std::unordered_set<long> nets_set(nets.begin(), nets.end());

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

    remove_free_nets(HC, 1);
    C.merge_free_vertices(HC);
    return HC;
}

} // namespace pmondriaan
