#pragma once

namespace pmondriaan {
	
/**
 * 
 */
class options {
	public:
		
		int sample_size() {return sample_size_; }
		
	private:
		int k_;
		int p_;
		
		double eps_;
		double eta_;
		
		std::string hypergraph_weigths_;
		std::string metric_;
		std::string bisection_mode_;
		
		int sample_size_;
		
		int coarsening_max_clustersize;
};
	
} // namespace pmondriaan