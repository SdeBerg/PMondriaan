#pragma once

#include <vector>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

namespace pmondriaan {

/**
 * List of vertices matched.
 */
class match_list {
	public:
		match_list()  { matches_ = std::vector<std::pair<int,int>>();}
		void add_match(int proc, int match) { matches_.push_back(std::make_pair(proc, match)); }
		
	
	private:
		std::vector<std::pair<int, int>> matches_;
};

/**
 * A contraction contains the information of how the hypergraph was contracted.
 */
class contraction {
	public:
		contraction(size_t size)  {
			matches_ = std::vector<pmondriaan::match_list>(size, pmondriaan::match_list());
		}
		
		void add_match(int sample, int match, int proc) { matches_[sample].add_match(proc, match); }
		
	
	private:
		std::vector<pmondriaan::match_list> matches_;
};
	



	
} // namespace pmondriaan