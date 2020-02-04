#include <iostream>

#include <pmondriaan.hpp>


int main () {	
	
	
	auto H = pmondriaan::read_hypergraph("../test/data/matrices/dolphins/dolphins.mtx", "degree");
	
	auto labels = pmondriaan::label_propagation(H, 5, 4);
	
	for (auto i = 0u; i < H.size(); i++) {
		std::cout << "i: " << i << "label: " << labels[i] << "\n";
	}
	
	return 0;
}