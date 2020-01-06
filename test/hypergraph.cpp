#include <iostream>
#include <fstream>

#include <pmondriaan.hpp>

int main () {	
	
	auto hypergraph = pmondriaan::read_hypergraph("../test/data/matrices/Trec5/Trec5.mtx", 3, 8);
	for (auto v : hypergraph) {
		std::cout << "id: " << v.id() << " weight: " << v.weight() << "\n nets: ";
		auto nets = v.nets();
		for (auto n : nets) {
			std::cout << n;
		}
		std::cout << "\n\n";
	}
	
	return 0;
}