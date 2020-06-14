#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <pmondriaan.hpp>


int main() {

    std::cout << "Creating a random hypergraph. Enter E V and L:\n";
    long E, V, L;
    std::cin >> E >> V >> L;
    pmondriaan::create_random_hypergraph("../test/data/matrices/random/nz" +
                                         std::to_string(L) + ".mtx",
                                         V, E, L);
    return 0;
}