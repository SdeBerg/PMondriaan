#pragma once

#include <limits>
#include <random>
#include <unordered_set>
#include <vector>

#include "hypergraph/hypergraph.hpp"
#include "util/interval.hpp"

namespace pmondriaan {

/**
 * The buckets containing the vertices of one part at the correct gain values.
 * Finding the vertex with highest gain can be done in O(1) time. Inserting and
 * removing vertices can also be done in O(1) time.
 */
class gain_buckets {
  public:
    gain_buckets(int size) {
        buckets = std::vector<std::unordered_set<int>>(size);
        max_index_present = size;
        min_value_present = std::numeric_limits<long>::max();
    }

    void insert(int id, long gain);

    // removes the element id from its bucket, returns true if this element existed and false otherwise
    bool remove(int id, long gain);

    int next();

    long gain_next();

    long gain_to_index(long gain) { return buckets.size() / 2 + gain; }
    long index_to_gain(long index) { return index - buckets.size() / 2; }

    void set_max_present(long index) { max_index_present = index; }

    void print();

  private:
    std::vector<std::unordered_set<int>> buckets;
    long max_index_present;
    long min_value_present;

    long find_next_index();
};

/**
 * This structure keeps track of the gain values during a KLFM pass. Each
 * unlocked vertex is contained in the buckets belonging to the part it is in.
 */
class gain_structure {
  public:
    gain_structure(pmondriaan::hypergraph& H, std::vector<std::vector<long>>& C)
    : H_(H), C_(C) {
        buckets =
        std::vector<pmondriaan::gain_buckets>(2, gain_buckets(compute_size_buckets()));
        gains = std::vector<long>(H.size());
        init_();
    }


    int part_next(long max_extra_weight_0, long max_extra_weight_1, std::mt19937& rng);

    int next(int part) { return buckets[part].next(); }

    long gain_next(int part) { return buckets[part].gain_next(); }

    void move(int v);

    void remove(int v);

    bool done();

    void add_gain(int v, long value);

  private:
    void init_();

    pmondriaan::hypergraph& H_;
    std::vector<std::vector<long>>& C_;
    std::vector<pmondriaan::gain_buckets> buckets;
    std::vector<long> gains;

    long compute_size_buckets();
};

} // namespace pmondriaan
