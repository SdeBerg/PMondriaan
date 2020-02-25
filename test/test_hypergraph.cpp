#include <iostream>
#include <fstream>
#include <array>
#include <vector>

#include <pmondriaan.hpp>


int main () {	
	
	
	auto hypergraph = pmondriaan::read_hypergraph("../test/data/matrices/cage3/cage3.mtx", "degree");
	
	hypergraph.print();
	
	return 0;
}