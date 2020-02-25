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
 * Uncoarsens the hypergraph HC into the hypergraph H.
 */
void uncoarsen_hypergraph(bulk::world& world,
                          pmondriaan::hypergraph& HC,
                          pmondriaan::hypergraph& H,
                          pmondriaan::contraction& C);

} // namespace pmondriaan
