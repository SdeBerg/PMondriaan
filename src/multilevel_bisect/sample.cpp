#include <vector>

#include "multilevel_bisect/sample.hpp"
#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

/**
 * Returns a vector of s randomly selected sample vertices. This is given as a list of indices.
 */
std::vector<int> sample_random(pmondriaan::hypergraph& H, int s) {
	auto samples = std::vector<int>();
	double size = (double)H.size();
	double s_left = (double)s;
	int current = 0;
	
	while(s_left > 0) {
		if (((double) rand() / (RAND_MAX)) < (s_left)/size) {
			s_left = s_left - 1.0;
			samples.push_back(current);
		}
		size = size - 1.0;
		current++;
	}
	
	return samples;
}

} // namespace pmondriaan
