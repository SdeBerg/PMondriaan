#include <cassert>
#include <vector>

#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/label_propagation.hpp"
#include "multilevel_bisect/sample.hpp"
#include "options.hpp"

namespace pmondriaan {

/**
 * Returns a vector of s randomly selected sample vertices. This is given as a list of indices.
 */
std::vector<long> sample_random(pmondriaan::hypergraph& H, long ns, std::mt19937& rng) {
    auto samples = std::vector<long>();
    double size = (double)H.size();

    assert(size > 0);
    assert(rng.max() > 0);

    double s_left = (double)ns;
    long current = 0;

    while (s_left > 0) {
        assert(size > 0.0);
        if (((double)rng() / rng.max()) < s_left / size) {
            s_left = s_left - 1.0;
            samples.push_back(current);
        }
        size = size - 1.0;
        current++;
    }

    return samples;
}


/**
 * Returns a vector of ns samples seleccted using the label propagation algorithm.
 */
std::vector<long>
sample_lp(pmondriaan::hypergraph& H, pmondriaan::options& opts, std::mt19937& rng) {

    auto labels =
    pmondriaan::label_propagation(H, opts.sample_size, opts.lp_max_iterations, 1, rng);
    auto count_label = std::vector<double>(opts.sample_size, 0.0);
    for (auto l : labels) {
        count_label[l]++;
    }

    long empty_count = 0;
    for (auto count : count_label) {
        if (count == 0) {
            empty_count++;
        }
    }

    auto samples = std::vector<long>();
    auto found_sample = std::vector<bool>(opts.sample_size, false);
    size_t number_samples_found = 0;
    auto current = 0u;
    while ((number_samples_found < opts.sample_size - empty_count) &&
           (current < H.size())) {
        long l = labels[current];
        if (!found_sample[l] &&
            (((double)rng() / (rng.max())) < 1.0 / (double)count_label[l])) {
            found_sample[l] = true;
            number_samples_found++;
            samples.push_back(current);
        }
        count_label[l]--;
        current++;
    }

    return samples;
}

} // namespace pmondriaan
