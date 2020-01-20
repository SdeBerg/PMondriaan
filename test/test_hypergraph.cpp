#include <iostream>
#include <fstream>
#include <array>

#include <pmondriaan.hpp>


int main () {	
	
	
	auto hypergraph = pmondriaan::read_hypergraph("../test/data/matrices/cage3/cage3.mtx");
	
	for (auto i = 0u; i < hypergraph.size(); i++) {
		auto v = hypergraph(i);
		std::cout << "id: " << v.id() << " weight: " << v.weight() << "\n";
		auto nets = v.nets();
		for (auto n : nets) {
			std::cout << n << "\n";
		}
	}
	
	return 0;
}