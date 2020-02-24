#include "hypergraph/readhypergraph.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <sstream>

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
pmondriaan::hypergraph read_hypergraph(std::string filename, std::string mode_weight) {
	
	int E, V;
	uint64_t L;
	
	std::ifstream fin(filename);
	if (fin.fail()){
		std::cerr << "Error: " << std::strerror(errno);
	}
	
	// Ignore headers and comments:
	while (fin.peek() == '%') {
		fin.ignore(2048, '\n');
	}
	
	// Read defining parameters:
	fin >> E >> V >> L;
	
	auto nets_list = std::vector<std::vector<int>>(V);
	auto vertex_list = std::vector<std::vector<int>>(E);
		
	// Read the data
	std::string line;
	while (std::getline(fin, line))
	{
		int e, v;
		std::istringstream iss(line);
		if (!(iss >> e >> v)) { break; }
		nets_list[v-1].push_back(e-1);
		vertex_list[e-1].push_back(v-1);
	}
	
	fin.close();
	
	auto vertices = std::vector<pmondriaan::vertex>();
	auto nets = std::vector<pmondriaan::net>();
  if (mode_weight == "one") {
    for (int i = 0; i < V; i++) {
		    vertices.push_back(pmondriaan::vertex(i, nets_list[i]));
		}
  }
	else if (mode_weight == "degree") {
    for (int i = 0; i < V; i++) {
		    vertices.push_back(pmondriaan::vertex(i, nets_list[i], nets_list[i].size()));
		}
	}

	for (int i = 0; i < E; i++) {
		nets.push_back(pmondriaan::net(i, vertex_list[i]));
	}
	
	return pmondriaan::hypergraph(V, vertices, nets);;
}

/**
 * Creates a distributed hypergraph from a graph in mtx format. 
 */
pmondriaan::hypergraph read_hypergraph(std::string filename, bulk::world& world, std::string mode_weight) {
	
	int s = world.rank();
    int p = world.active_processors();
	
	auto E = bulk::var<int>(world);
	auto V = bulk::var<int>(world);
	auto L = bulk::var<uint64_t>(world);
	
	
	/* Queue used to assign vertices to a processor. A message is given
	   as a vertex id and a net id that contains the vertex. */
	auto vertices_queue = bulk::queue<int, int>(world);
	
	if (s == 0) {
		
		std::ifstream fin(filename);
		if (fin.fail()){
			std::cerr << "Error: " << std::strerror(errno);
		}
		
		// Ignore headers and comments:
		while (fin.peek() == '%') {
			fin.ignore(2048, '\n');
		}
		
		// Read defining parameters:
		fin >> E >> V >> L;
	
		E.broadcast(E);
		V.broadcast(V);
		L.broadcast(L);
		
		fin.close();
	}
	
	world.sync();
	
	auto partitioning = bulk::block_partitioning<1>({V}, {p});
	
	if (s == 0) {
		std::ifstream fin(filename);
		if (fin.fail()){
			std::cerr << "Error: " << std::strerror(errno);
		}
		
		while (fin.peek() == '%') {
			fin.ignore(2048, '\n');
		}
		
		// Read the data
		std::string line;
		while (std::getline(fin, line))
		{
			int e, v;
			std::istringstream iss(line);
			if (!(iss >> e >> v)) { break; } // error
			vertices_queue(partitioning.owner({v-1})).send(v-1, e-1);
		}
		fin.close();
	}
	
	world.sync();
	
	// List of nets for each vertex
	auto nets_list = std::vector<std::vector<int>>(partitioning.local_count(s));
	auto vertex_list = std::vector<std::vector<int>>(E);
	for (const auto& [v,n] : vertices_queue) {
		int v_loc = partitioning.local({v})[0];
		nets_list[v_loc].push_back(n);
		vertex_list[n].push_back(v);
	}
	
	auto vertices = std::vector<pmondriaan::vertex>();
	auto nets = std::vector<pmondriaan::net>();
  if (mode_weight == "one") {
    for (int i = 0; i < partitioning.local_count(s); i++) {
			vertices.push_back(pmondriaan::vertex(partitioning.global({i}, s)[0], nets_list[i]));
		}
  }
	else if (mode_weight == "degree") {
    for (int i = 0; i < partitioning.local_count(s); i++) {
			vertices.push_back(pmondriaan::vertex(partitioning.global({i}, s)[0], nets_list[i], nets_list[i].size()));
		}
  }
  else {
    std::cerr << "Error: unknown mode_weight";
  }
	for (int i = 0; i < E; i++) {
		nets.push_back(pmondriaan::net(i, vertex_list[i]));
	}
	
	auto H = pmondriaan::hypergraph(V, vertices, nets);
	
	pmondriaan::remove_free_nets(world, H);
	
	return H;
}

} // namespace pmondriaan
