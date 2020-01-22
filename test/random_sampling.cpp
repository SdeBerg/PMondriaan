#include <iostream>
#include <fstream>
#include <vector>

#include <pmondriaan.hpp>


int main () {	
	
	
	auto hypergraph = pmondriaan::read_hypergraph("../test/data/matrices/gemat11/gemat11.mtx", "degree");
	
	auto samples = pmondriaan::sample_random(hypergraph, 10);
	
	for (auto s : samples) {
		std::cout << s << "\n";
	}
	
	return 0;
}