#include <vector>

#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/label_propagation.hpp"
#include "multilevel_bisect/sample.hpp"
#include "options.hpp"

namespace pmondriaan {

/**
 * Returns a vector of s randomly selected sample vertices. This is given as a list of indices.
 */
std::vector<int> sample_random(pmondriaan::hypergraph& H, int ns) {
    auto samples = std::vector<int>();
    double size = (double)H.size();
    double s_left = (double)ns;
    int current = 0;

    while (s_left > 0) {
        if (((double)rand() / (RAND_MAX)) < (s_left) / size) {
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
std::vector<int> sample_lp(pmondriaan::hypergraph& H, pmondriaan::options& opts) {

    auto labels = pmondriaan::label_propagation(H, opts.sample_size, opts.lp_max_iterations, 1);
    auto count_label = std::vector<double>(opts.sample_size, 0.0);
    for (auto l : labels) {
        count_label[l]++;
    }

    int empty_count = 0;
    for (auto count : count_label) {
        if (count == 0) {
            empty_count++;
        }
    }

    auto samples = std::vector<int>();
    auto found_sample = std::vector<bool>(opts.sample_size, false);
    int number_samples_found = 0;
    auto current = 0u;
    while (number_samples_found < opts.sample_size - empty_count && current < H.size()) {
        int l = labels[current];
        if (!found_sample[l] &&
            ((double)rand() / (RAND_MAX)) < 1.0 / (double)count_label[l]) {
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
