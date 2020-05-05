#include <algorithm>
#include <cmath>
#include <random>
#include <stack>
#include <stdlib.h>
#include <string>
#include <unordered_set>
#include <vector>

#include <iostream>

#include "algorithm.hpp"
#include "bisect.hpp"
#include "hypergraph/hypergraph.hpp"
#include "options.hpp"
#include "recursive_bisection.hpp"
#include "util/interval.hpp"
#include "work_item.hpp"

namespace pmondriaan {

/**
 * Recursively bisects a hypergraph into k parts with imbalance parameter epsilon.
 * Eta is the load imbalance parameter for the imbalance over the processors during computation.
 */
void recursive_bisect(bulk::world& world,
                      pmondriaan::hypergraph& H,
                      long k,
                      double epsilon,
                      double eta,
                      pmondriaan::options& opts) {

    auto s = world.rank();
    auto p = world.active_processors();

    std::mt19937 rng(s + 1);

    auto global_weight = pmondriaan::global_weight(world, H);
    long maxweight = ((1.0 + epsilon) * global_weight) / k;

    auto sub_world = world.split(0);

    // start and end index of the current part to be split
    long start = 0;
    long end = H.size();
    auto jobs = std::stack<pmondriaan::work_item>();

    interval labels = {0, k - 1};

    long weight_mypart = global_weight;

    if (p == 1) {
        jobs.push(pmondriaan::work_item(start, end, labels.low, labels.high, weight_mypart));
    }

    // while we need to give more than one label and need to use more than one processor, we bisect the hypergraph in parallel
    while ((labels.length() > 0) && (p > 1)) {

        long k_ = labels.length() + 1;

        long k_low = k_ / 2;
        long k_high = k_ - k_low;

        auto max_global_weights =
        compute_max_global_weight(k_, k_low, k_high, weight_mypart, maxweight);

        // part 0 will always have the smallest weight
        auto weight_parts = bisect(*sub_world, H, opts, max_global_weights[0],
                                   max_global_weights[1], start, end, labels, rng);

        auto total_weight_0 = bulk::sum(*sub_world, weight_parts[0]);
        auto total_weight_1 = bulk::sum(*sub_world, weight_parts[1]);

        if (sub_world->rank() == 0) {
            world.log("s: %d, weight part %d: %d, weight part %d: %d", world.rank(),
                      labels.low, total_weight_0, labels.high, total_weight_1);
        }

        // number of processors working on the low and high part respectively
        long p_low =
        (double)total_weight_0 / (double)(total_weight_0 + total_weight_1) * (double)p + 0.5;
        if ((k_low == 1) & (k_high > 1)) {
            p_low = 0;
        }
        long p_high = p - p_low;

        long my_part = 1;
        if (s < p_low) {
            my_part = 0;
        }

        // personal low and high label for the next round
        interval new_labels = {labels.low + my_part * k_low,
                               labels.high - (1 - my_part) * k_high};

        if (new_labels.length() > 0) {

            long new_max_local_weight;
            if (my_part == 0) {
                new_max_local_weight =
                ceil((double)total_weight_0 / (double)p_low) * (1.0 + eta);
            } else {
                new_max_local_weight =
                ceil((double)total_weight_1 / (double)p_high) * (1.0 + eta);
            }

            /**
             * Redistribute the hypergraph over the processors such that all vertices with label_low are on
             * processors with my_part 0 and with label_high on processors with my_part 1.
             */
            end = redistribute_hypergraph(*sub_world, H, my_part, labels.low,
                                          labels.high, new_max_local_weight,
                                          weight_parts[0], weight_parts[1], p_low);

            // the processors in the new world will work together on the next part
            sub_world = sub_world->split(my_part);
            p = sub_world->active_processors();
            s = sub_world->rank();

            if (my_part == 0) {
                weight_mypart = total_weight_0;
            } else {
                weight_mypart = total_weight_1;
            }

            // if we are the only processor left, we add a sequential job
            if (p == 1) {
                jobs.push(pmondriaan::work_item(start, end, new_labels.low,
                                                new_labels.high, weight_mypart));
            }
        }

        sub_world->sync();

        labels = new_labels;
    }

    // we do the rest of the work sequentially
    while (!jobs.empty()) {
        auto job = jobs.top();
        jobs.pop();

        start = job.start();
        end = job.end();
        labels = {job.label_low(), job.label_high()};
        weight_mypart = job.weight();

        while (labels.length() > 0) {
            long k_ = labels.length() + 1;

            long k_low = k_ / 2;
            long k_high = k_ - k_low;

            auto max_global_weights =
            compute_max_global_weight(k_, k_low, k_high, weight_mypart, maxweight);

            auto weight_parts = bisect(*sub_world, H, opts, max_global_weights[0],
                                       max_global_weights[1], start, end, labels, rng);

            interval labels_0 = {labels.low, labels.high - k_high};
            interval labels_1 = {labels.low + k_low, labels.high};

            long end_1 = end;

            if (labels_1.length() > 0) {
                reorder_hypergraph(H, start, end, labels.low, labels.high);
                jobs.push(pmondriaan::work_item(end, end_1, labels_1.low,
                                                labels_1.high, weight_parts[1]));
            }
            weight_mypart = weight_parts[0];
            labels = labels_0;
        }
    }

    world.sync();
}

std::vector<long>
compute_max_global_weight(long k_, long k_low, long k_high, long weight_mypart, long maxweight) {

    double eps = (double)(maxweight * k_) / (double)weight_mypart - 1.0;

    long q_low = std::log2(k_low) + 1;
    long q_high = std::log2(k_high) + 1;

    double delta_low = eps / q_low;
    double delta_high = eps / q_high;

    long global_weight_low =
    ((double)weight_mypart / (double)k_) * (double)k_low * (1.0 + delta_low);
    long global_weight_high =
    ((double)weight_mypart / (double)k_) * (double)k_high * (1.0 + delta_high);

    return {global_weight_low, global_weight_high};
}

/**
 * Redistributes the hypergraph such that processors with my_part 0 contain all
 * vertices with label_low and all with my_part 1 contain all vertices with
 * label_high. Returns the end of part 0 if part 1 is not assigned any
 * processors.
 */
long redistribute_hypergraph(bulk::world& world,
                             pmondriaan::hypergraph& H,
                             long my_part,
                             long label_low,
                             long label_high,
                             long max_local_weight,
                             long weight_part_0,
                             long weight_part_1,
                             long p_low) {

    long surplus_0 = (1 - my_part) * max_local_weight - weight_part_0;
    long surplus_1 = my_part * max_local_weight - weight_part_1;
    auto all_surplus_0 = bulk::gather_all(world, surplus_0);
    auto all_surplus_1 = bulk::gather_all(world, surplus_1);

    /* we move all vertices of part 0 to a separate vector, be be put at the
       back of the vertices later if part 0 is finished */
    auto vertices_0 = std::vector<pmondriaan::vertex>();

    // queue for all received vertices
    auto q = bulk::queue<long, long, long, long[]>(world);
    // the weights of new nets are communicated through this queue
    auto cost_queue = bulk::queue<long, long>(world);
    reduce_surplus(world, H, label_high, all_surplus_1, q, cost_queue);

    if (p_low > 0) {
        reduce_surplus(world, H, label_low, all_surplus_0, q, cost_queue);
    } else {
        long i = H.size();
        while (i >= 0) {
            if (H(i).part() == label_low) {
                std::move(H.vertices().begin() + i, H.vertices().begin() + i + 1,
                          std::back_inserter(vertices_0));
                H.vertices().erase(H.vertices().begin() + i);
            }
            i--;
        }
    }

    world.sync();

    // We add the new nets we received
    for (const auto& [net_id, cost] : cost_queue) {
        H.add_net(net_id, std::vector<long>(), cost);
    }

    for (const auto& [id, weight, part, nets] : q) {
        H.vertices().push_back({id, nets, weight});
        H.vertices().back().set_part(part);
        H.add_to_nets(H.vertices().back());
    }

    if (p_low == 0) {
        H.vertices().insert(end(H.vertices()), begin(vertices_0), end(vertices_0));
    }

    H.update_map();

    return H.size() - vertices_0.size();
}

/**
 * Distributes the surplus vertices of a processor over the other processors.
 * Returns 0 if all surplus is gone, -1 if there is still surplus left.
 */
long reduce_surplus(bulk::world& world,
                    pmondriaan::hypergraph& H,
                    long label,
                    bulk::coarray<long>& surplus,
                    bulk::queue<long, long, long, long[]>& q,
                    bulk::queue<long, long>& cost_queue) {

    auto s = world.rank();
    auto p = world.active_processors();

    if (surplus[s] >= 0) {
        return 0;
    }

    // index of next vertex to be sent
    long index = 0;
    long surplus_others = 0;
    long t = s + 1;
    while (surplus[s] < 0) {
        if (t >= p) {
            t = 0;
        }
        if (t == s) {
            std::cerr << "Error: failed to lose all surplus\n";
            return -1;
        }
        surplus_others += surplus[t];
        if (surplus[t] > 0 && surplus_others > 0) {

            auto nets_to_send = std::unordered_set<long>();

            long max = std::min(surplus[t], surplus_others);
            max = std::min(max, -1 * surplus[s]);

            long total_sent = 0;
            while (total_sent < max) {
                // search for the next vertex to be sent
                while (H(index).part() != label) {
                    index++;
                }

                auto v = H(index);

                total_sent += v.weight();
                nets_to_send.insert(v.nets().begin(), v.nets().end());
                q(t).send(v.id(), v.weight(), v.part(), v.nets());

                long id = v.id();
                H.remove_from_nets(id);
                H.map().erase(id);
                H.vertices().erase(H.vertices().begin() + index);
            }

            // We send the information on the cost of the nets.
            for (auto n : nets_to_send) {
                cost_queue(t).send(n, H.net(n).cost());
            }
            surplus[s] += total_sent;
            surplus_others = 0;
        }
        t++;
    }

    return 0;
}

/**
 * Reorders the hypergraph such that all vertices with label_high are at the end of the vertex list.
 */
void reorder_hypergraph(pmondriaan::hypergraph& H, long start, long& end, long label_low, long label_high) {
    long pivot = start;
    while (pivot < end) {
        if (H(pivot).part() == label_high) {
            while ((H(end - 1).part() != label_low) && (end - 1 != pivot)) {
                end--;
            }
            std::swap(H(pivot), H(end - 1));
        }
        pivot++;
    }
    end--;
    H.update_map();
}

} // namespace pmondriaan
