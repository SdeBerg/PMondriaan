#include "hypergraph/readhypergraph.hpp"

#include <cstring>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

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
std::optional<pmondriaan::hypergraph>
read_hypergraph_istream(std::istream& fin, std::string mode_weight) {

    int E, V;
    uint64_t L;

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
    std::getline(fin, line);
    while (std::getline(fin, line)) {
        int e, v;
        std::istringstream iss(line);
        if (!(iss >> e >> v)) {
            return std::nullopt;
        }
        nets_list[v - 1].push_back(e - 1);
        vertex_list[e - 1].push_back(v - 1);
    }

    // fin.close();

    auto vertices = std::vector<pmondriaan::vertex>();
    auto nets = std::vector<pmondriaan::net>();
    if (mode_weight == "one") {
        for (int i = 0; i < V; i++) {
            vertices.push_back(pmondriaan::vertex(i, nets_list[i]));
        }
    } else if (mode_weight == "degree") {
        for (int i = 0; i < V; i++) {
            vertices.push_back(pmondriaan::vertex(i, nets_list[i], nets_list[i].size()));
        }
    }

    for (int i = 0; i < E; i++) {
        if (!vertex_list[i].empty()) {
            nets.push_back(pmondriaan::net(i, vertex_list[i]));
        }
    }

    auto H = pmondriaan::hypergraph(V, E, vertices, nets);
    remove_free_nets(H);
    return std::move(H);
}

/**
 * Creates a distributed hypergraph from a graph in mtx format.
 */
std::optional<pmondriaan::hypergraph>
read_hypergraph_istream(std::istream& fin, bulk::world& world, std::string mode_weight) {

    int s = world.rank();
    int p = world.active_processors();

    int E, V;
    uint64_t L;

    // Ignore headers and comments:
    while (fin.peek() == '%') {
        fin.ignore(2048, '\n');
    }

    // Read defining parameters:
    fin >> E >> V >> L;

    auto partitioning = bulk::block_partitioning<1>({V}, {p});

    // List of nets for each vertex
    auto nets_list = std::vector<std::vector<int>>(partitioning.local_count(s));
    auto vertex_list = std::vector<std::vector<int>>(E);

    // Read the data
    std::string line;
    std::getline(fin, line);
    for (auto i = 0u; i < L; i++) {
        std::getline(fin, line);
        int e, v;
        std::istringstream iss(line);
        if (!(iss >> e >> v)) {
            return std::nullopt;
        } // error
        if (partitioning.owner({v - 1}) == s) {
            int v_loc = partitioning.local({v - 1})[0];
            nets_list[v_loc].push_back(e - 1);
            vertex_list[e - 1].push_back(v - 1);
        }
    }
    world.sync();

    auto vertices = std::vector<pmondriaan::vertex>();
    auto nets = std::vector<pmondriaan::net>();
    if (mode_weight == "one") {
        for (int i = 0; i < partitioning.local_count(s); i++) {
            vertices.push_back(
            pmondriaan::vertex(partitioning.global({i}, s)[0], nets_list[i]));
        }
    } else if (mode_weight == "degree") {
        for (int i = 0; i < partitioning.local_count(s); i++) {
            vertices.push_back(pmondriaan::vertex(partitioning.global({i}, s)[0],
                                                  nets_list[i], nets_list[i].size()));
        }
    } else {
        std::cerr << "Error: unknown mode_weight";
    }
    for (int i = 0; i < E; i++) {
        if (!vertex_list[i].empty()) {
            nets.push_back(pmondriaan::net(i, vertex_list[i]));
        }
    }

    auto H = pmondriaan::hypergraph(V, E, vertices, nets);
    pmondriaan::remove_free_nets(world, H);
    world.log("global0 %d", H.global_number_nets());
    world.sync();
    return std::move(H);
}

std::optional<pmondriaan::hypergraph> read_hypergraph(std::string file, std::string mode_weight) {
    std::ifstream fs(file);
    if (fs.fail()) {
        std::cerr << "Error: " << std::strerror(errno);
        return std::nullopt;
    }
    return read_hypergraph_istream(fs, mode_weight);
}

std::optional<pmondriaan::hypergraph>
read_hypergraph(std::string file, bulk::world& world, std::string mode_weight) {
    std::ifstream fs(file);
    if (fs.fail()) {
        std::cerr << "Error: " << std::strerror(errno);
        return std::nullopt;
    }
    return read_hypergraph_istream(fs, world, mode_weight);
}

} // namespace pmondriaan
