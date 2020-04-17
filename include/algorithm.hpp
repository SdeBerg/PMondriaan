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
 * Find owner of minimal value of x
 */
template <typename T>
int owner_min(bulk::var<T>& x) {
    auto& world = x.world();

    auto images = bulk::gather_all(world, x.value());

    T smallest = images[0];
    int best_proc = 0;
    for (int t = 1; t < world.active_processors(); ++t) {
        if (images[t] <= smallest) {
            smallest = images[t];
            best_proc = t;
        }
    }
    return best_proc;
}

} // namespace pmondriaan