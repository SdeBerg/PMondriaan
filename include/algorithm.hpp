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
 * Perform a left-associative fold over a coarray.
 */
template <typename T, typename Func, typename S = T>
std::vector<S> foldl(bulk::coarray<T>& x, Func f, S start_value = {}) {
    auto& world = x.world();
    auto result = std::vector<T>(x.size(), start_value);
	
	for (auto i = 0u; i < x.size(); i++) {
		auto images = bulk::gather_all(world, x[i]);
		
		for (int t = 0; t < world.active_processors(); ++t) {
			// apply f iteratively to the current value, and each remote value
			f(result[i], images[t]);
		}
	}
	
    return result;
}
	
} // namespace pmondriaan