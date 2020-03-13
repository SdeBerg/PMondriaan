#include <fstream>
#include <iostream>
#include <vector>

#include <pmondriaan.hpp>


int main() {

    auto hypergraph =
    pmondriaan::read_hypergraph("../test/data/matrices/dolphins/dolphins.mtx", "one");
    if (!hypergraph) {
        std::cerr << "Error: failed to load hypergraph\n";
        return -1;
    }
    auto H = hypergraph.value();

    std::mt19937 rng(127);
    auto L = label_propagation_bisect(H, 100, 40, 40, rng);
    for (auto l : L) {
        std::cout << l << "\n";
    }
    return 0;
}