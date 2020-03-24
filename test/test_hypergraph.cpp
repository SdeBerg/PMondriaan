#include <array>
#include <fstream>
#include <iostream>
#include <vector>

#include <pmondriaan.hpp>


int main() {

    auto H =
    pmondriaan::read_hypergraph("../test/data/matrices/dolphins/dolphins.mtx", "one");
    auto hypergraph = H.value();
    hypergraph.print();
    std::random_device random;
    std::mt19937 rng(random());
    auto L = pmondriaan::label_propagation_bisect(hypergraph, 100, 36, 36, rng);

    for (auto i = 0u; i < hypergraph.size(); i++) {
        std::cout << "i: " << i << " label: " << L[i] << "\n";
    }

    return 0;
}