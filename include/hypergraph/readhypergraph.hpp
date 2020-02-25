#pragma once

#include <cstring>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include <hypergraph/hypergraph.hpp>

namespace pmondriaan {

/**
 * Creates a hypergraph from a graph in mtx format.
 */
pmondriaan::hypergraph read_hypergraph(std::string filename, std::string mode_weight = "one");

/**
 * Creates a distributed hypergraph from a graph in mtx format. The mode states
 * how the weights of the vertices are computed.
 */
pmondriaan::hypergraph
read_hypergraph(std::string filename, bulk::world& world, std::string mode_weight = "one");

} // namespace pmondriaan
