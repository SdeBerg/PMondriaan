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
pmondriaan::hypergraph read_hypergraph(std::string filename);

/**
 * Creates a distributed hypergraph from a graph in mtx format. 
 */
pmondriaan::hypergraph read_hypergraph(std::string filename, bulk::world& world);

} // namespace pmondriaan
