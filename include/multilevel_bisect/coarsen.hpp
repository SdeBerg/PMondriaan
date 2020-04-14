#include <random>
#include <string>
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

namespace pmondriaan {

/**
 * Coarsens the hypergraph H and returns a hypergraph HC in parallel.
 */
pmondriaan::hypergraph coarsen_hypergraph_par(bulk::world& world,
                                              pmondriaan::hypergraph& H,
                                              pmondriaan::contraction& C,
                                              pmondriaan::options& opts,
                                              std::mt19937& rng);


/**
 * Sends match request to the owners of the best matches found using the
 * improduct computation. Returns the local matches.
 */
void request_matches(pmondriaan::hypergraph& H,
                     pmondriaan::contraction& C,
                     bulk::queue<int, long, int[]>& sample_queue,
                     bulk::queue<int, int>& accepted_matches,
                     const std::vector<int>& indices_samples,
                     pmondriaan::options& opts);

/**
 * First merges the nets and weight of all vertices matched to a sample and
 * then sends this information to the owner of the sample.
 */
void send_information_matches(pmondriaan::hypergraph& H,
                              bulk::queue<int, int>& accepted_matches,
                              bulk::queue<int, long, int[], long[]>& info_queue,
                              std::vector<bool>& matched,
                              int sample_size);

pmondriaan::hypergraph contract_hypergraph(bulk::world& world,
                                           pmondriaan::hypergraph& H,
                                           const std::vector<int> samples,
                                           bulk::queue<int, long, int[], long[]>& matches,
                                           std::vector<bool>& matched);

/**
 * Coarsens the hypergraph H and returns a hypergraph HC sequentially.
 */
pmondriaan::hypergraph coarsen_hypergraph_seq(bulk::world& world,
                                              pmondriaan::hypergraph& H,
                                              pmondriaan::contraction& C,
                                              pmondriaan::options& opts,
                                              std::mt19937& rng);

/**
 * Add a copy of a vertex v to a list of vertices.
 */
void add_v_to_list(std::vector<pmondriaan::vertex>& v_list, pmondriaan::vertex& v);

pmondriaan::hypergraph contract_hypergraph(bulk::world& world,
                                           pmondriaan::hypergraph& H,
                                           pmondriaan::contraction& C,
                                           std::vector<std::vector<int>>& matches,
                                           std::vector<pmondriaan::vertex>& new_vertices);

} // namespace pmondriaan
