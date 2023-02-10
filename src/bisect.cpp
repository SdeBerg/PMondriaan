#include <limits>
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
constexpr long stopping_time_par = 3;
}
constexpr bool print_time = true;
constexpr bool simplify_duplicates = true;

namespace pmondriaan {

/**
 * Bisect a hypergraph using the given bisection method and returns the weights of the two parts.
 */
std::vector<long> bisect(bulk::world& world,
                         pmondriaan::hypergraph& H,
                         pmondriaan::options& opts,
                         long max_weight_0,
                         long max_weight_1,
                         long start,
                         long end,
                         interval labels,
                         std::mt19937& rng) {

    auto weight_parts = std::vector<long>(2);
    auto p = world.active_processors();
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
                                long start,
                                long end,
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
                                    long start,
                                    long end,
                                    interval labels,
                                    std::mt19937& rng) {

    // a hypergraph containing only vertices with indices between start and end is created
    auto H_reduced = pmondriaan::create_new_hypergraph(world, H, start, end);

    // the number of parallel coursenings performed
    size_t nc_par = 0;

    auto HC_list = std::vector<pmondriaan::hypergraph>{H_reduced};
    auto C_list = std::vector<pmondriaan::contraction>();
    C_list.push_back({});

    if (world.active_processors() > 1) {
        C_list[0].merge_free_vertices(world, HC_list[0]);
    } else {
        C_list[0].merge_free_vertices(HC_list[0]);
    }

    auto time = bulk::util::timer();

    // PARALLEL COARSENING PHASE
    if (world.active_processors() > 1) {
        // size of the hypergraph at which we stop the pararallel coarsening process
        size_t coarsening_nrvertices_par =
        std::max(opts.coarsening_nrvertices, world.active_processors() * opts.sample_size *
                                             parameters::stopping_time_par);

        double ratio = 1.0;
        while ((HC_list[nc_par].global_size() > coarsening_nrvertices_par) &&
               (nc_par < opts.coarsening_maxrounds) && ratio > 0.05) {

            C_list.push_back({});
            HC_list.push_back(coarsen_hypergraph_par(world, HC_list[nc_par],
                                                     C_list[nc_par + 1], opts, rng));

            nc_par++;
            ratio =
            (double)(HC_list[nc_par - 1].global_size() - HC_list[nc_par].global_size()) /
            (double)HC_list[nc_par - 1].global_size();

            if (world.rank() == 0) {
                if (print_time) {
                    world.log("s: %d, time in iteration par coarsening: %lf",
                              world.rank(), time.get_change());
                }
                world.log("After iteration %d, size is %d (par)", nc_par,
                          HC_list[nc_par].global_size());
            }
        }

        // we now communicate the entire hypergraph to all processors using a queue containing id, weight and nets
        auto vertex_queue = bulk::queue<long, long, long[]>(world);
        for (auto& v : HC_list[nc_par].vertices()) {
            for (long t = 0; t < world.rank(); t++) {
                vertex_queue(t).send(v.id(), v.weight(), v.nets());
            }
            for (long t = world.rank() + 1; t < world.active_processors(); t++) {
                vertex_queue(t).send(v.id(), v.weight(), v.nets());
            }
        }
        // we also communicate the cost of all local nets
        auto net_cost_queue = bulk::queue<long, long>(world);
        for (auto& n : HC_list[nc_par].nets()) {
            for (long t = 0; t < world.rank(); t++) {
                net_cost_queue(t).send(n.id(), n.cost());
            }
            for (long t = world.rank() + 1; t < world.active_processors(); t++) {
                net_cost_queue(t).send(n.id(), n.cost());
            }
        }
        HC_list.push_back(HC_list[nc_par]);
        C_list.push_back({});
        nc_par++;
        world.sync();

        if (print_time && (world.rank() == 0)) {
            world.log("s: %d, time in sending vertices and net costs: %lf",
                      world.rank(), time.get_change());
        }

        // the received nets and vertices are added to the coarsened hypergraph
        for (const auto& [net_id, cost] : net_cost_queue) {
            HC_list[nc_par].add_net(net_id, std::vector<long>(), cost);
        }
        for (const auto& [id, weight, nets] : vertex_queue) {
            HC_list[nc_par].add_vertex(id, nets, weight);
            HC_list[nc_par].add_to_nets(HC_list[nc_par].vertices().back());
        }

        if (print_time && (world.rank() == 0)) {
            world.log("s: %d, time in creating new hypergraph: %lf",
                      world.rank(), time.get_change());
        }
    }

    time.get();

    auto nc_tot = nc_par;
    auto max_rounds = opts.coarsening_maxrounds;
    if (world.active_processors() > 1) {
        max_rounds++;
    }

    // SEQUENTIAL COARSENING PHASE
    while ((HC_list[nc_tot].global_size() > opts.coarsening_nrvertices) &&
           (nc_tot < max_rounds)) {
        C_list.push_back({});
        time.get();

        if (simplify_duplicates) {
            simplify_duplicate_nets(HC_list[nc_tot]);
            if (world.rank() == 0) {
                if (print_time) {
                    world.log("s: %d, time seq simplifying duplicate nets: %lf",
                              world.rank(), time.get_change());
                    time.get();
                }
            }
        }

        HC_list.push_back(coarsen_hypergraph_seq(world, HC_list[nc_tot],
                                                 C_list[nc_tot + 1], opts, rng));

        nc_tot++;
        if (world.rank() == 0) {
            if (print_time) {
                world.log("s: %d, time in iteration seq coarsening: %lf",
                          world.rank(), time.get_change());
            }
            world.log("After iteration %d, size is %d (seq)", nc_tot - 1,
                      HC_list[nc_tot].global_size());
        }
    }

    time.get();
    // INITIAL PARTITIONING PHASE
    auto cut = pmondriaan::initial_partitioning(HC_list[nc_tot], max_weight_0,
                                                max_weight_1, opts, rng);

    if (world.rank() == 0) {
        if (print_time) {
            world.log("s: %d, time in initial partitioning: %lf", world.rank(),
                      time.get_change());
        }
        world.log("s %d: cut after initial partitioning: %d", world.rank(), cut);
    }

    // SEQUENTIAL UNCOARSENING PHASE
    while (nc_tot > nc_par) {
        nc_tot--;
        cut = pmondriaan::uncoarsen_hypergraph_seq(HC_list[nc_tot + 1], HC_list[nc_tot],
                                                   C_list[nc_tot + 1], opts, max_weight_0,
                                                   max_weight_1, cut, rng);

        HC_list.pop_back();
        C_list.pop_back();

        if (world.rank() == 0) {
            if (print_time) {
                world.log("s: %d, time in iteration seq uncoarsening: %lf",
                          world.rank(), time.get_change());
            }
            world.log("s %d: cut after seq uncoarsening: %d", world.rank(), cut);
        }

        if (simplify_duplicates) {
            time.get();
            HC_list[nc_tot].reset_duplicate_nets();
            if (world.rank() == 0) {
                if (print_time) {
                    world.log("s: %d, time in seq resetting duplicate nets: "
                              "%lf",
                              world.rank(), time.get_change());
                }
            }
        }
    }

    // PARALLEL UNCOARSENING PHASE
    if (world.active_processors() > 1) {
        // we find the best solution so far of all partitioners
        bulk::var<long> cut_size(world);
        cut_size = cut;
        if (HC_list[nc_par].weight_part(0) > max_weight_0 ||
            HC_list[nc_par].weight_part(1) > max_weight_1) {
            cut_size = std::numeric_limits<long>::max();
        }
        auto best_proc = pmondriaan::owner_min(cut_size);

        // the processor that has found the best solution now sends the labels to all others
        auto label_queue = bulk::queue<long, long>(world);
        if (world.rank() == best_proc) {
            cut_size.broadcast(cut);
            for (auto& v : HC_list[nc_par].vertices()) {
                for (long t = 0; t < world.active_processors(); t++) {
                    label_queue(t).send(v.id(), v.part());
                }
            }
        }
        world.sync();

        // the local coarsened hypergraph is given the received labels
        nc_par--;
        for (const auto& [id, part] : label_queue) {
            if (HC_list[nc_par].is_local(id)) {
                HC_list[nc_par](HC_list[nc_par].local_id(id)).set_part(part);
            }
        }

        cut = cut_size;
        if (world.rank() == 0) {
            world.log("s %d: cut at start par uncoarsening: %d", world.rank(), cut);
        }
        time.get();
        while (nc_par > 0) {
            nc_par--;
            cut = pmondriaan::uncoarsen_hypergraph_par(world, HC_list[nc_par + 1],
                                                       HC_list[nc_par],
                                                       C_list[nc_par + 1], opts, max_weight_0,
                                                       max_weight_1, cut, rng);

            HC_list.pop_back();
            C_list.pop_back();

            if (world.rank() == 0) {
                if (print_time) {
                    world.log("s: %d, time in iteration par uncoarsening: %lf",
                              world.rank(), time.get_change());
                }
                world.log("s %d: cut after par uncoarsening: %d", world.rank(), cut);
            }
        }
    }

    if (world.active_processors() > 1) {
        C_list[0].assign_free_vertices(world, HC_list[0], max_weight_0, max_weight_1, rng);
    } else {
        C_list[0].assign_free_vertices(HC_list[0], max_weight_0, max_weight_1, rng);
    }

    // the vertices of the original hypergraph are labelled using the labeling of the reduced hypergraph
    for (auto& v : HC_list[0].vertices()) {
        H(H.local_id(v.id())).set_part(labels(v.part()));
    }

    auto weight_parts = std::vector<long>(2);
    weight_parts[0] = HC_list[0].weight_part(0);
    weight_parts[1] = HC_list[0].weight_part(1);
    return weight_parts;
}

} // namespace pmondriaan
