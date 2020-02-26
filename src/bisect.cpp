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

namespace pmondriaan {

/**
 * Bisect a hypergraph using the given bisection method and returns the weights of the two parts.
 */
std::vector<long> bisect(bulk::world& world,
                         pmondriaan::hypergraph& H,
                         std::string bisect_mode,
                         std::string sampling_mode,
                         pmondriaan::options& opts,
                         std::string metric,
                         long max_weight_0,
                         long max_weight_1,
                         int start,
                         int end,
                         interval labels) {

    auto weight_parts = std::vector<long>(2);
    int p = world.active_processors();
    if (bisect_mode == "random") {
        weight_parts = bisect_random(world, H, max_weight_0 / p, max_weight_1 / p,
                                     start, end, labels);
    }

    if (bisect_mode == "multilevel") {
        weight_parts = bisect_multilevel(world, H, opts, sampling_mode, metric, max_weight_0,
                                         max_weight_1, start, end, labels);
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
                                interval labels) {

    std::mt19937 rng(world.rank() + 1);

    auto max_weight_parts = std::vector<long>{max_weight_0, max_weight_1};

    auto weight_parts = std::vector<long>(2);

    for (auto i = start; i < end; i++) {
        int part = rng() % 2;
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
                                    interval labels) {

    auto H_reduced = pmondriaan::create_new_hypergraph(world, H, start, end);

    long nc_par = 0;

    auto HC_list = std::vector<pmondriaan::hypergraph>();
    auto C_list = std::vector<pmondriaan::contraction>();
    HC_list.push_back(H_reduced);

    while ((HC_list[nc_par].global_size() > opts.coarsening_nrvertices) &&
           (nc_par < opts.coarsening_maxrounds)) {
        C_list.push_back(pmondriaan::contraction());
        HC_list.push_back(coarsen_hypergraph(world, HC_list[nc_par], C_list[nc_par], opts, sampling_mode));
        nc_par++;
        world.log("After iteration %d, size is %d", nc_par, HC_list[nc_par].global_size());
	}
  
	
	//TODO: communicate hypergraph to all processors and coarsen further sequentially
	/*long nc_tot = nc_par + 1;
	while ((HC_list[nc_tot].global_size() > opts.coarsening_nrvertices) &&
           (nc_tot < opts.coarsening_maxrounds)) {}*/

    pmondriaan::initial_partitioning(world, HC_list[nc_par], max_weight_0,
                                     max_weight_1, labels);

	//TODO: uncoarsen locally best solution sequentially up to the parallel point, communicate best solution, pick best solution and receive info on this
    
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
