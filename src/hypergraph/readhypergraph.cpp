#include "hypergraph/readhypergraph.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
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
pmondriaan::hypergraph read_hypergraph(std::string filename) {
	
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
	for (auto l = 0u; l < L; l++) {
		int e, v;
		double data;
		fin >> e >> v >> data;
		nets_list[v-1].push_back(e-1);
		vertex_list[e-1].push_back(v-1);
	}
	
	fin.close();
	
	auto vertices = std::vector<pmondriaan::vertex>();
	auto nets = std::vector<pmondriaan::net>();
	for (int i = 0; i < V; i++) {
		vertices.push_back(pmondriaan::vertex(i, nets_list[i]));
	}
	for (int i = 0; i < E; i++) {
		nets.push_back(pmondriaan::net(i, vertex_list[i]));
	}
	
	return pmondriaan::hypergraph(V, vertices, nets);;
}

/**
 * Creates a distributed hypergraph from a graph in mtx format. 
 */
pmondriaan::hypergraph read_hypergraph(std::string filename, bulk::world& world) {
	
	int s = world.rank();
    int p = world.active_processors();
	
	auto E = bulk::var<int>(world);
	auto V = bulk::var<int>(world);
	auto L = bulk::var<uint64_t>(world);
	
	
	// Queue where vertices will be assigned to a processor. A message is given
	// as a vertex id and a net id that contains the vertex.
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
		
		// Reading the first line again:
		fin >> E >> V >> L;
		
		// Read the data
		for (auto l = 0u; l < L; l++) {
			int e, v;
			double data;
			fin >> e >> v >> data;
			vertices_queue(partitioning.owner({v-1})).send(v-1, e-1);
		}
		fin.close();
	}
	
	world.sync();
	
	// List of nets for each vertex
	auto nets_list = std::vector<std::vector<int>>(partitioning.local_count(s));
	auto vertex_list = std::vector<std::vector<int>>(E);
	for (auto vertex_net : vertices_queue) {
		int v, n;
		std::tie(v, n) = vertex_net;
		int v_loc = partitioning.local({v})[0];
		nets_list[v_loc].push_back(n);
		vertex_list[n].push_back(v_loc);
	}
	
	auto vertices = std::vector<pmondriaan::vertex>();
	auto nets = std::vector<pmondriaan::net>();
	for (int i = 0; i < partitioning.local_count(s); i++) {
		vertices.push_back(pmondriaan::vertex(partitioning.global({i}, s)[0], nets_list[i]));
	}
	for (int i = 0; i < E; i++) {
		nets.push_back(pmondriaan::net(i, vertex_list[i]));
	}
	
	return pmondriaan::hypergraph(V, vertices, nets);
}

} // namespace pmondriaan