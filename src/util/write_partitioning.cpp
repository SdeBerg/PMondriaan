#include <fstream>
#include <iostream>
#include <string>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

bool partitioning_to_file(bulk::world& world, pmondriaan::hypergraph& H, std::string file) {
    if (world.rank() == 0) {
        std::ofstream out(file);
        if (out.fail()) {
            std::cerr << "Error: " << std::strerror(errno);
            return false;
        }
        out << H.global_size() << "\n";
        out.close();
    }
    for (int turn = 0; turn < world.active_processors(); turn++) {
        if (turn == world.rank()) {
            std::ofstream out(file, std::ios::app);
            if (out.fail()) {
                std::cerr << "Error: " << std::strerror(errno);
                return false;
            }
            for (auto& v : H.vertices()) {
                out << v.id() + 1 << " " << v.part() + 1 << "\n";
            }
            out.close();
        }
        world.sync();
    }
    return true;
}

} // namespace pmondriaan