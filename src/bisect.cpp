#include <math.h>
#include <random>
#include <stdlib.h>
#include <vector>

#include "bisect.hpp"
#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/coarsen.hpp"
#include "multilevel_bisect/initial_partitioning.hpp"
#include "multilevel_bisect/uncoarsen.hpp"
#include "options.hpp"
#include "util/interval.hpp"
#include "algorithm.hpp"

namespace pmondriaan {

/**
 * Bisect a hypergraph using the given bisection method and returns the weights of the two parts.
 */
std::vector<long> bisect(bulk::world& world,pmondriaan::hypergraph& H,std::string bisect_mode,std::string sampling_mode,
                         pmondriaan::options& opts,
                         std::string metric,
                         long max_weight_0,
                         long max_weight_1,
                         int start,
                         int end,
                         interval labels,
                         std::mt19937& rng) {

    auto weight_parts = std::vector<long>(2);
    int p = world.active_processors();
    if (bisect_mode == "random") {
        weight_parts = bisect_random(world, H, max_weight_0 / p, max_weight_1 / p,
                                     start, end, labels, rng);
    }

    if (bisect_mode == "multilevel") {
        weight_parts = bisect_multilevel(world, H, opts, sampling_mode, metric, max_weight_0,
                                         max_weight_1, start, end, labels, rng);
    }

    return weight_parts;
}


/**
 * Randomly bisects a hypergraph under the balance constraint and returns the weights of the two parts.
 */
std::vector<long> bisect_random(bulk::world& world,
                                pmondriaan::hypergraph& H,
                                long max_weight_0,
                                long max_weight_1,
                                int start,
                                int end,
                                interval labels,
                                std::mt19937& rng) {

    auto max_weight_parts = std::vector<long>{max_weight_0, max_weight_1};

    auto weight_parts = std::vector<long>(2);
    for (auto i = start; i < end; i++) {
        int part = rng() % 2;
		//world.log("weight %d, sum %d, max %d, part %d", H(i).weight(), weight_parts[part] + H(i).weight(), max_weight_parts[part], part);
        // if the max weight is exceeded, we add the vertex to the other part
        if (weight_parts[part] + H(i).weight() > max_weight_parts[part]) {
            part = (part + 1) % 2;
        }

        if (part == 0) {
            H(i).set_part(labels.low);
        } else {
            H(i).set_part(labels.high);
        }
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
                                    std::string sampling_mode,
                                    std::string metric,
                                    long max_weight_0,
                                    long max_weight_1,
                                    int start,
                                    int end,
                                    interval labels,
                                    std::mt19937& rng) {

    auto H_reduced = pmondriaan::create_new_hypergraph(world, H, start, end);

    long nc_par = 0;

    auto HC_list = std::vector<pmondriaan::hypergraph>();
    auto C_list = std::vector<pmondriaan::contraction>();
    HC_list.push_back(H_reduced);

	auto coarsening_nrvertices_par = std::max(opts.coarsening_nrvertices, world.active_processors() * opts.sample_size * 5);

	while ((HC_list[nc_par].global_size() > coarsening_nrvertices_par) &&
           (nc_par < opts.coarsening_maxrounds)) {
        C_list.push_back(pmondriaan::contraction());
        HC_list.push_back(coarsen_hypergraph_par(world, HC_list[nc_par], C_list[nc_par], opts, sampling_mode, rng));
        nc_par++;
        world.log("After iteration %d, size is %d (par)", nc_par, HC_list[nc_par].global_size());
	}

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

	world.sync();

	for (const auto& [id, weight, nets] : vertex_queue) {
		HC_list[nc_par].vertices().push_back({id, nets, weight});
        HC_list[nc_par].add_to_nets(HC_list[nc_par].vertices().back());
	}
	HC_list[nc_par].update_map();

	long nc_tot = nc_par;
	while ((HC_list[nc_tot].global_size() > opts.coarsening_nrvertices) &&
           (nc_tot < opts.coarsening_maxrounds)) {
        C_list.push_back(pmondriaan::contraction());
        HC_list.push_back(coarsen_hypergraph_seq(world, HC_list[nc_tot], C_list[nc_tot], opts, rng));
        nc_tot++;
        world.log("After iteration %d, size is %d (seq)", nc_tot, HC_list[nc_tot].global_size());
	}

    pmondriaan::initial_partitioning(world, HC_list[nc_tot], max_weight_0,
                                     max_weight_1, labels, rng);

	while (nc_tot > nc_par) {
        nc_tot--;
        pmondriaan::uncoarsen_hypergraph(world, HC_list[nc_tot + 1], HC_list[nc_tot], C_list[nc_tot]);
    }

	//we find the best solution of all partitioners
	bulk::var<long> cut(world);
	cut = pmondriaan::cutsize(HC_list[nc_par], metric);
	auto best_proc = pmondriaan::smallest(cut);

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
        pmondriaan::uncoarsen_hypergraph(world, HC_list[nc_par + 1], HC_list[nc_par], C_list[nc_par]);
    }

    for (auto& v : HC_list[0].vertices()) {
        H(H.local_id(v.id())).set_part(v.part());
    }

    auto weight_parts = std::vector<long>(2);
    weight_parts[0] = HC_list[0].weight_part(labels.low);
    weight_parts[1] = HC_list[0].weight_part(labels.high);

    return weight_parts;
}

} // namespace pmondriaan
