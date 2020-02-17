#include <vector>
#include <string>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "bisect.hpp"
#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/sample.hpp"
#include "hypergraph/contraction.hpp"

namespace pmondriaan {

/**
 * Coarses the hypergraph H and returns a hypergraph HC.
 */
pmondriaan::hypergraph coarsen_hypergraph(bulk::world& world, pmondriaan::hypergraph& H, pmondriaan::options& opts, std::string sampling_mode);


/** Sends match request to the owners of the best matches found using the improduct computation.
 * Returns the local matches. 
 */
void request_matches(pmondriaan::hypergraph& H, auto& sample_queue, auto& accepted_matches, const std::vector<int>& indices_samples, pmondriaan::options& opts);

pmondriaan::hypergraph contract_hypergraph(bulk::world& world, pmondriaan::hypergraph& H, const std::vector<int> samples, auto& matches, std::vector<bool>& matched);


} // namespace pmondriaan
