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
 * Create a coarray with images holding the given value on each processor in the given range.
 *
 * This function takes an argument, and writes it to the appropriate element on
 * each remote processor.
 */
template <typename T>
bulk::coarray<T> gather_all(bulk::world& world, std::vector<int> range, T value) {
    bulk::coarray<T> xs(world, range[1] - range[0]);

    for (int t = range[0]; t < range[1]; ++t) {
        xs(t)[world.rank() - range[0]] = value;
    }
	
    world.sync();
	
    return xs;
}

/**
 * Perform a left-associative fold for a range of processors
 */
template <typename T, typename Func, typename S = T>
S foldl(bulk::var<T>& x, Func f, std::vector<int> range, S start_value = {}) {
    auto& world = x.world();
    auto result = start_value;
	
    auto images = pmondriaan::gather_all(world, range, x.value());

    for (int t = 0; t < range[1] - range[0]; ++t) {
        // apply f iteratively to the current value, and each remote value
        f(result, images[t]);
    }
    return result;
}
	
/**
 * Perform a left-associative fold over a coarray.
 */
template <typename T, typename Func, typename S = T>
std::vector<T> foldl(bulk::coarray<T>& x, Func f, S start_value = {}) {
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