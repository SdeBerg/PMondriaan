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
                     bulk::queue<long, long, long[]>& sample_queue,
                     bulk::queue<long, long>& accepted_matches,
                     const std::vector<long>& indices_samples,
                     pmondriaan::options& opts);

/**
 * First merges the nets and weight of all vertices matched to a sample and
 * then sends this information to the owner of the sample.
 */
void send_information_matches(bulk::world& world,
                              pmondriaan::hypergraph& H,
                              bulk::queue<long, long>& accepted_matches,
                              bulk::queue<long, long, long[], long[]>& info_queue,
                              std::vector<bool>& matched,
                              long sample_size);

pmondriaan::hypergraph contract_hypergraph(bulk::world& world,
                                           pmondriaan::hypergraph& H,
                                           pmondriaan::contraction& C,
                                           const std::vector<long> samples,
                                           bulk::queue<long, long, long[], long[]>& matches,
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
                                           std::vector<std::vector<long>>& matches,
                                           std::vector<pmondriaan::vertex>& new_vertices);

} // namespace pmondriaan
