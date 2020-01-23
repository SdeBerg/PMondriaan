#include <iostream>
#include <fstream>
#include <vector>

#include <pmondriaan.hpp>


int main () {	
	
	
	auto hypergraph = pmondriaan::read_hypergraph("../test/data/matrices/gemat11/gemat11.mtx", "degree");
	
	pmondriaan::bisect_multilevel(world, hypergraph, 0, 0, 0, hypergraph.size());
	
	return 0;
}