#include <vector>
#include <algorithm>
#include <iostream>

#include "multilevel_bisect/label_propagation.hpp"
#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

/**
 * Performs label propagation to create l groups of labels on the hypergraph H.
 */
std::vector<int> label_propagation(pmondriaan::hypergraph& H, int l, int max_iter) {
	//counts of all labels for each net
	auto C = std::vector<std::vector<long>>(H.nets().size(), std::vector<long>(l, 0));
	//the labels of the vertices
	auto L = std::vector<int>(H.size());
	
	for (auto& label : L) {
		int random = rand();
		label = random % l;
	}
	
	std::vector<int> indices(H.size()); 
	std::iota(indices.begin(), indices.end(), 0); 
	
	auto T = std::vector<long>(l);
	bool change = true;
	int iterations = 0;
	
	while (change && (iterations < max_iter)) {
		std::random_shuffle(indices.begin(), indices.end());
		change = false;
		for (auto i : indices) {
			std::fill(T.begin(), T.end(), 0);
			
			/* First compute the sum of the counts */
			for (auto n : H(i).nets()) {
				C[n][L[i]]--;
				for (int j = 0; j < l; j++) {
					T[j] += C[n][j];
				}
			}
			
			/* Now compute argmax(T), where we break ties randomly */
			long max = -1;
			auto label_max = std::vector<int>();
			for (int j = 0; j < l; j++) {
				if (T[j] >= max) {
					max = T[j];
					label_max.push_back(j);
				}
			}
			
			/* Set new label and change if it is changed */
			int random = rand();
			int new_label = label_max[random % label_max.size()];
			if (new_label != L[i]) {
				change = true;
				L[i] = new_label;
			}
			
			/* Update the counts */
			for (auto n : H(i).nets()) {
				C[n][L[i]]++;
			}
		}
		iterations++;
	}
	
	return L;
}

} // namespace pmondriaan