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
 * Perform a left-associative fold over a coarray for each element.
 *
 * This function applies a function to the images of the elements of a coarray.
 *
 * The cost is `F * p * X + p * g * X + l`, where `F` is the number of
 * flops performed during a single call to `f` and `X` is the size of the coarray.
 *
 * \tparam T the type of the value held by \c x.
 * \tparam Func the binary function to apply to the images of \c x.
 *
 * \param x the coarray to fold over
 * \param f a binary function that takes two arguments of type `T`.
 *
 * \returns a vector with the result of the expression
 *              \f[ f(f(f(f(x(0)[i], x(1)[i]), x(2)[i]), ...), x(p-1)[i]). \f]
 *          for each i, which is computed at each core.
 */
template <typename T, typename Func, typename S = T>
std::vector<T> foldl_each(bulk::coarray<T>& x, Func f, S start_value = {}) {
    auto& world = x.world();
    auto result = std::vector<T>(x.size(), start_value);

    bulk::coarray<T> images(world, world.active_processors() * x.size());

    for (int t = 0; t < world.active_processors(); ++t) {
        for (auto i = 0u; i < x.size(); ++i) {
            images(t)[i + world.rank() * x.size()] = x[i];
        }
    }

    world.sync();

    for (auto i = 0u; i < x.size(); i++) {
        for (int t = 0; t < world.active_processors(); ++t) {
            // apply f iteratively to the current value, and each remote value
            f(result[i], images[t * x.size() + i]);
        }
    }

    return result;
}


/**
 * Perform a left-associative fold over a coarray.
 *
 * This function applies a function to the images of a coarray.
 *
 * The cost is `F * (p + X) + p * g + l`, where `F` is the number of
 * flops performed during a single call to `f` and `X` the size of the coarray.
 *
 * \tparam T the type of the value held by \c x.
 * \tparam Func the binary function to apply to the images of \c x.
 *
 * \param x the coarray to fold over
 * \param f a binary function that takes two arguments of type `T`.
 *
 * \returns the result of the expression
 *              \f[ f(f(f(f(x(0)[o], x(0)[1]), x(0)[2]), ...), x(p-1)[X - 1]).
 * \f] which is computed at each core.
 */
template <typename T, typename Func, typename S = T>
S foldl(bulk::coarray<T>& x, Func f, S start_value = {}) {
    auto& world = x.world();
    auto local_result = start_value;
    auto result = start_value;

    for (auto i = 0u; i < x.size(); ++i) {
        // apply f iteratively to each local value
        f(local_result, x[i]);
    }

    auto images = bulk::gather_all(world, local_result);

    for (int t = 0; t < world.active_processors(); ++t) {
        // apply f iteratively to the current value, and each remote value
        f(result, images[t]);
    }
    return result;
}

/**
 * Find owner of smallest value of x
 */
template <typename T>
int smallest(bulk::var<T>& x) {
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