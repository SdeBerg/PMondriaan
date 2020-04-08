#pragma once

#include <optional>
#include <string>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include <hypergraph/hypergraph.hpp>
#include <optional>

namespace pmondriaan {

/**
 * Creates a hypergraph from a matrix in mtx format.
 */
std::optional<pmondriaan::hypergraph>
read_hypergraph_istream(std::istream& fin, std::string mode_weight = "one");

/**
 * Creates a distributed hypergraph from a graph in mtx format. The mode states
 * how the weights of the vertices are computed.
 */
std::optional<pmondriaan::hypergraph>
read_hypergraph_istream(std::istream& fin, bulk::world& world, std::string mode_weight = "one");

/**
 * Creates a hypergraph from a file that contains a matrix in mtx format.
 */
std::optional<pmondriaan::hypergraph>
read_hypergraph(std::string file, std::string mode_weight = "one");

/**
 * Creates a distributed hypergraph from a file that contains a matrix in mtx format.
 */
std::optional<pmondriaan::hypergraph>
read_hypergraph(std::string file, bulk::world& world, std::string mode_weight = "one");

} // namespace pmondriaan
