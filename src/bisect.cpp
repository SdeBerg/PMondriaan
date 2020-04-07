#include <math.h>
#include <random>
#include <stdlib.h>
#include <vector>

#include "algorithm.hpp"
#include "bisect.hpp"
#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/coarsen.hpp"
#include "multilevel_bisect/initial_partitioning.hpp"
#include "multilevel_bisect/uncoarsen.hpp"
#include "options.hpp"
#include "util/interval.hpp"

namespace parameters {
constexpr int stopping_time_par = 5;
}

namespace pmondriaan {

/**
 * Bisect a hypergraph using the given bisection method and returns the weights of the two parts.
 */
std::vector<long> bisect(bulk::world& world,
                         pmondriaan::hypergraph& H,
                         pmondriaan::options& opts,
                         long max_weight_0,
                         long max_weight_1,
                         int start,
                         int end,
                         interval labels,
                         std::mt19937& rng) {

    auto weight_parts = std::vector<long>(2);
    int p = world.active_processors();
    if (opts.bisection_mode == pmondriaan::bisection::random) {
        weight_parts =
        bisect_random(H, max_weight_0 / p, max_weight_1 / p, start, end, labels, rng);
    }

    if (opts.bisection_mode == pmondriaan::bisection::multilevel) {
        weight_parts = bisect_multilevel(world, H, opts, max_weight_0,
                                         max_weight_1, start, end, labels, rng);
    }

    return weight_parts;
}


/**
 * Randomly bisects a hypergraph under the balance constraint and returns the weights of the two parts.
 */
std::vector<long> bisect_random(pmondriaan::hypergraph& H,
                                long max_weight_0,
                                long max_weight_1,
                                int start,
                                int end,
                                interval labels,
                                std::mt19937& rng) {

    auto max_weight_parts = std::vector<long>{max_weight_0, max_weight_1};

    auto weight_parts = std::vector<long>(2);
    for (auto i = start; i < end; i++) {
        auto part = rng() % 2;
        // if the max weight is exceeded, we add the vertex to the other part
        if (weight_parts[part] + H(i).weight() > max_weight_parts[part]) {
            part = (part + 1) % 2;
        }

        H(i).set_part(labels(part));
        weight_parts[part] += H(i).weight();
    }

    return weight_parts;
}

/**
 * Bisects a hypergraph using the multilevel framework.
 */
std::vector<long> bisect_multilevel(bulk::world& world,
                                    pmondriaan::hypergraph& H,
                                    pmondriaan::options& opts,
                                    long max_weight_0,
                                    long max_weight_1,
                                    int start,
                                    int end,
                                    interval labels,
                                    std::mt19937& rng) {
    auto H_reduced = pmondriaan::create_new_hypergraph(world, H, start, end);

    long nc_par = 0;

    auto HC_list = std::vector<pmondriaan::hypergraph>{H_reduced};
    auto C_list = std::vector<pmondriaan::contraction>();

    auto time = bulk::util::timer();

    if (world.active_processors() > 1) {
        auto coarsening_nrvertices_par =
        std::max(opts.coarsening_nrvertices, world.active_processors() * opts.sample_size *
                                             parameters::stopping_time_par);

        while ((HC_list[nc_par].global_size() > coarsening_nrvertices_par) &&
               (nc_par < opts.coarsening_maxrounds)) {
            C_list.push_back({});
            HC_list.push_back(coarsen_hypergraph_par(world, HC_list[nc_par],
                                                     C_list[nc_par], opts, rng));
            nc_par++;
            world.log("After iteration %d, size is %d (par)", nc_par,
                      HC_list[nc_par].global_size());
        }

        world.log("s: %d, time in par coarsening: %lf", world.rank(), time.get_change());
        // we now communicate the entire hypergraph to all processors using a queue containing id, weight and nets
        auto vertex_queue = bulk::queue<int, long, int[]>(world);
        for (auto& v : HC_list[nc_par].vertices()) {
            for (int t = 0; t < world.rank(); t++) {
                vertex_queue(t).send(v.id(), v.weight(), v.nets());
            }
            for (int t = world.rank() + 1; t < world.active_processors(); t++) {
                vertex_queue(t).send(v.id(), v.weight(), v.nets());
            }
        }
        world.log("s: %d, time in sending vertices: %lf", world.rank(), time.get_change());
        HC_list.push_back(HC_list[nc_par]);
        C_list.push_back({});
        nc_par++;
        world.sync();

        for (const auto& [id, weight, nets] : vertex_queue) {
            HC_list[nc_par].vertices().push_back({id, nets, weight});
            HC_list[nc_par].add_to_nets(HC_list[nc_par].vertices().back());
        }
        world.log("s: %d, time in creating new hypergraph: %lf", world.rank(),
                  time.get_change());
        // TODO: Take this out in the final version! This is only included to get reproducibility
        sort(HC_list[nc_par].vertices().begin(), HC_list[nc_par].vertices().end(),
             [](pmondriaan::vertex a, pmondriaan::vertex b) {
                 return a.id() < b.id();
             });
        world.log("s: %d, time in sorting: %lf", world.rank(), time.get_change());
        HC_list[nc_par].update_map();
    }
    time.get();
    long nc_tot = nc_par;
    while ((HC_list[nc_tot].global_size() > opts.coarsening_nrvertices) &&
           (nc_tot - 1 < opts.coarsening_maxrounds)) {
        C_list.push_back(pmondriaan::contraction());
        time.get();
        HC_list.push_back(coarsen_hypergraph_seq(world, HC_list[nc_tot],
                                                 C_list[nc_tot], opts, rng));
        world.log("s: %d, time in iteration seq coarsening: %lf", world.rank(),
                  time.get_change());
        nc_tot++;
        world.log("After iteration %d, size is %d (seq)", nc_tot - 1,
                  HC_list[nc_tot].global_size());
    }
    // world.log("s: %d, time in sequential coarsening: %lf", world.rank(), time.get_change());
    auto cut = pmondriaan::initial_partitioning(HC_list[nc_tot], max_weight_0,
                                                max_weight_1, opts, rng);

    while (nc_tot > nc_par) {
        nc_tot--;
        cut = pmondriaan::uncoarsen_hypergraph_seq(world, HC_list[nc_tot + 1],
                                                   HC_list[nc_tot], C_list[nc_tot],
                                                   opts, max_weight_0,
                                                   max_weight_1, cut, rng);
    }

    if (world.active_processors() > 1) {
        // we find the best solution of all partitioners
        bulk::var<long> cut_size(world);
        cut_size = cut;
        auto best_proc = pmondriaan::owner_min(cut_size);

        // the processor that has found the best solution now sends the labels to all others
        auto label_queue = bulk::queue<int, int>(world);
        if (world.rank() == best_proc) {
            for (auto& v : HC_list[nc_par].vertices()) {
                for (int t = 0; t < world.active_processors(); t++) {
                    label_queue(t).send(v.id(), v.part());
                }
            }
        }
        world.sync();

        if (world.rank() != best_proc) {
            for (const auto& [id, part] : label_queue) {
                HC_list[nc_par](HC_list[nc_par].local_id(id)).set_part(part);
            }
        }

        while (nc_par > 0) {
            nc_par--;
            pmondriaan::uncoarsen_hypergraph(world, HC_list[nc_par + 1],
                                             HC_list[nc_par], C_list[nc_par]);
        }
    }

    for (auto& v : HC_list[0].vertices()) {
        H(H.local_id(v.id())).set_part(labels(v.part()));
    }

    auto weight_parts = std::vector<long>(2);
    weight_parts[0] = HC_list[0].weight_part(0);
    weight_parts[1] = HC_list[0].weight_part(1);

    return weight_parts;
}

} // namespace pmondriaan
