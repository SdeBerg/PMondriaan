#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"
#include "util/write_partitioning.hpp"

namespace pmondriaan {

bool partitioning_to_file(bulk::world& world, pmondriaan::hypergraph& H, std::string file, int k) {
    auto starts = start_parts(world, H, k);
    if (world.rank() == 0) {
        std::ofstream out(file);
        if (out.fail()) {
            std::cerr << "Error: " << std::strerror(errno);
            return false;
        }
        out << "%%MatrixMarket distributed-matrix coordinate pattern general\n";
        out << H.global_number_nets() << " " << H.global_size() << " "
            << H.nr_nz() << " " << k << "\n";
        for (auto i : starts) {
            out << i << "\n";
        }
        out.close();
    }
    for (int i = 0; i < k; i++) {
        for (long turn = 0; turn < world.active_processors(); turn++) {
            if (turn == world.rank()) {
                std::ofstream out(file, std::ios::app);
                if (out.fail()) {
                    std::cerr << "Error: " << std::strerror(errno);
                    return false;
                }
                for (auto& v : H.vertices()) {
                    if (v.part() == i) {
                        for (auto n : v.nets()) {
                            out << n + 1 << " " << v.id() + 1 << "\n";
                        }
                    }
                }
                out.close();
            }
            world.sync();
        }
    }
    return true;
}

std::vector<long> start_parts(bulk::world& world, pmondriaan::hypergraph& H, int k) {
    auto counts = bulk::coarray<long>(world, k);
    for (int i = 0; i < k; i++) {
        counts[i] = 0;
    }
    for (auto& v : H.vertices()) {
        counts[v.part()] += v.degree();
    }
    auto result =
    bulk::foldl_each(counts, [](auto& lhs, auto rhs) { lhs += rhs; });
    result.insert(result.begin(), 0);
    for (int i = 1; i < k + 1; i++) {
        result[i + 1] += result[i];
    }
    return result;
}

} // namespace pmondriaan