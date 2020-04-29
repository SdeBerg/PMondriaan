#include <random>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/contraction.hpp"
#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/sample.hpp"

namespace pmondriaan {

/**
 * Uncoarsens the hypergraph HC sequentially into the hypergraph H.
 * The cutsize is then optimized using the KLFM algorithm. Returns
 * the cutsize of the partitioning found.
 */
long uncoarsen_hypergraph_seq(pmondriaan::hypergraph& HC,
                              pmondriaan::hypergraph& H,
                              pmondriaan::contraction& C,
                              pmondriaan::options& opts,
                              long max_weight_0,
                              long max_weight_1,
                              long cut_size,
                              std::mt19937& rng);

/**
 * Uncoarsens the hypergraph HC into the hypergraph H.
 * The cutsize is then optimized using the parallel KLFM algorithm. Returns
 * the cutsize of the partitioning found.
 */
long uncoarsen_hypergraph_par(bulk::world& world,
                              pmondriaan::hypergraph& HC,
                              pmondriaan::hypergraph& H,
                              pmondriaan::contraction& C,
                              pmondriaan::options& opts,
                              long max_weight_0,
                              long max_weight_1,
                              long cut_size,
                              std::mt19937& rng);

/**
 * Uncoarsens the hypergraph HC into the hypergraph H.
 */
void uncoarsen_hypergraph(bulk::world& world,
                          pmondriaan::hypergraph& HC,
                          pmondriaan::hypergraph& H,
                          pmondriaan::contraction& C);

/**
 * Uncoarsens the hypergraph HC into the hypergraph H.
 */
void uncoarsen_hypergraph(pmondriaan::hypergraph& HC,
                          pmondriaan::hypergraph& H,
                          pmondriaan::contraction& C);

} // namespace pmondriaan
