#include <algorithm>
#include <iostream>
#include <vector>

#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/label_propagation.hpp"

namespace pmondriaan {

/**
 * Performs label propagation to create l groups of labels on the hypergraph H.
 */
std::vector<int>
label_propagation(pmondriaan::hypergraph& H, int l, int max_iter, int min_size, std::mt19937& rng) {
    // counts of all labels for each net
    auto C = std::vector<std::vector<long>>(H.nets().size(), std::vector<long>(l, 0));
    auto size_L = std::vector<long>(l, 0);
    // the labels of the vertices
    auto L = std::vector<int>(H.size());

    for (auto i = 0u; i < L.size(); i++) {
        auto random = rng();
        L[i] = random % l;
        size_L[L[i]]++;
        for (auto n : H(i).nets()) {
            C[H.local_id_net(n)][L[i]]++;
        }
    }

    std::vector<int> indices(H.size());
    std::iota(indices.begin(), indices.end(), 0);

    auto T = std::vector<long>(l);
    bool change = true;
    int iterations = 0;

    while (change && (iterations < max_iter)) {
        std::shuffle(indices.begin(), indices.end(), rng);
        change = false;
        for (auto i : indices) {
            std::fill(T.begin(), T.end(), 0);

            if (size_L[L[i]] > min_size) {
                // First compute the sum of the counts
                for (auto n : H(i).nets()) {
                    C[H.local_id_net(n)][L[i]]--;
                    for (int j = 0; j < l; j++) {
                        T[j] += C[H.local_id_net(n)][j];
                    }
                }

                // Now compute argmax(T), where we break ties randomly
                long max = -1;
                auto label_max = std::vector<int>();
                for (int j = 0; j < l; j++) {
                    if (T[j] == max) {
                        label_max.push_back(j);
                    } else if (T[j] > max) {
                        max = T[j];
                        label_max.clear();
                        label_max.push_back(j);
                    }
                }

                // Set new label and change if it is changed
                auto random = rng();
                int new_label = label_max[random % label_max.size()];
                if (new_label != L[i]) {
                    change = true;
                    size_L[L[i]]--;
                    L[i] = new_label;
                    size_L[L[i]]++;
                }

                // Update the counts
                for (auto n : H(i).nets()) {
                    C[H.local_id_net(n)][L[i]]++;
                }
            }
        }
        iterations++;
    }

    return L;
}

std::vector<int> label_propagation_bisect(pmondriaan::hypergraph& H,
                                          std::vector<std::vector<long>>& C,
                                          int max_iter,
                                          long max_weight_0,
                                          long max_weight_1,
                                          std::mt19937& rng) {

    auto L = std::vector<int>(H.size());
    // stores that weight that can still be assigned to the labels
    auto weight_L = std::vector<long>(2, 0);
    weight_L[0] = max_weight_0;
    weight_L[1] = max_weight_1;
    // the labels of the vertices

    for (auto i = 0u; i < H.size(); i++) {
        auto random = rng();
        L[i] = random % 2;
        if ((weight_L[L[i]] - H(i).weight()) < 0) {
            L[i] = (L[i] + 1) % 2;
        }
        weight_L[L[i]] -= H(i).weight();
        for (auto n : H(i).nets()) {
            C[H.local_id_net(n)][L[i]]++;
        }
    }

    std::vector<int> indices(H.size());
    std::iota(indices.begin(), indices.end(), 0);

    auto T = std::vector<long>(2);
    bool change = true;
    int iterations = 0;
    while (change && (iterations < max_iter)) {
        std::shuffle(indices.begin(), indices.end(), rng);
        change = false;
        for (auto i : indices) {
            std::fill(T.begin(), T.end(), 0);

            // First compute the sum of the counts
            for (auto n : H(i).nets()) {
                C[H.local_id_net(n)][L[i]]--;
                for (int j = 0; j < 2; j++) {
                    T[j] += C[H.local_id_net(n)][j];
                }
            }

            // Now compute argmax(T), where we break ties randomly
            long max = -1;
            auto label_max = std::vector<int>();
            for (int j = 0; j < 2; j++) {
                if (T[j] == max) {
                    label_max.push_back(j);
                } else if (T[j] > max) {
                    max = T[j];
                    label_max.clear();
                    label_max.push_back(j);
                }
            }

            // Set new label and change if it is changed
            if (max != -1) {
                auto random = rng();
                int new_label = label_max[random % label_max.size()];
                if ((new_label != L[i]) && ((weight_L[new_label] - H(i).weight()) >= 0)) {
                    change = true;
                    weight_L[L[i]] += H(i).weight();
                    L[i] = new_label;
                    weight_L[L[i]] -= H(i).weight();
                }
            }

            // Update the counts
            for (auto n : H(i).nets()) {
                C[H.local_id_net(n)][L[i]]++;
            }
        }
        iterations++;
    }

    return L;
}

} // namespace pmondriaan