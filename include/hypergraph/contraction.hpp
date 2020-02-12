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
 * A contraction contains the information of how the hypergraph was contracted.
 */
class contraction {
	public:
		contraction(size_t size)  {
			vertices_ = std::vector<std::vector<int>>(size);
		}
		
		
		
	
	private:
		std::vector<std::vector<int>> vertices_;
};
	
	
} // namespace pmondriaan