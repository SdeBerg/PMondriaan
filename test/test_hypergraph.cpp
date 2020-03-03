#include <iostream>
#include <fstream>
#include <array>
#include <vector>

#include <pmondriaan.hpp>


int main () {


	auto H = pmondriaan::read_hypergraph("../test/data/matrices/dolphins/dolphins.mtx", "one");
	auto hypergraph = H.value();
	hypergraph.print();

	auto L = pmondriaan::label_propagation_bisect(hypergraph, 100, 36, 36);

	for (auto i = 0u; i < hypergraph.size(); i++) {
		std::cout << "i: " << i << " label: " << L[i] <<"\n";
	}

	return 0;
}