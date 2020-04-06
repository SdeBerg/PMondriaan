#pragma once

#include <string>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

bool partitioning_to_file(bulk::world& world, pmondriaan::hypergraph& H, std::string file);

} // namespace pmondriaan