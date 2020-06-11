#include <fstream>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include "hypergraph/hypergraph.hpp"
#include "util/random_hypergraph.hpp"

namespace pmondriaan {

bool create_random_hypergraph(std::string file, long n, long m, long nz) {
    std::random_device rd;
    std::mt19937 rng(rd());

    std::ofstream out(file);
    if (out.fail()) {
        std::cerr << "Error: " << std::strerror(errno);
        return false;
    }
    out << "%%MatrixMarket matrix coordinate pattern general\n";
    out << m << " " << n << " " << nz << "\n";
    long nz_added = 0;
    auto list_j = std::vector<std::unordered_set<long>>(m);
    while (nz_added < nz) {
        auto i = rng() % m;
        auto j = rng() % n;
        if (list_j[i].insert(j).second) {
            out << i + 1 << " " << j + 1 << "\n";
            nz_added++;
        }
    }
    out.close();
    return true;
}

} // namespace pmondriaan