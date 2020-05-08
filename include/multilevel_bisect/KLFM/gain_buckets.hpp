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
    gain_buckets(long size) {
        buckets = std::vector<std::unordered_set<long>>(size);
        max_index_present = size;
        min_value_present = index_to_gain(0);
    }

    void insert(long id, long gain);

    // removes the element id from its bucket, returns true if this element existed and false otherwise
    bool remove(long id, long gain);

    long next();

    long gain_next();

    long gain_to_index(long gain) { return buckets.size() / 2 + gain; }
    long index_to_gain(long index) { return index - buckets.size() / 2; }

    void set_max_present(long index) { max_index_present = index; }

    void print();

  private:
    std::vector<std::unordered_set<long>> buckets;
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


    long part_next(long max_extra_weight_0, long max_extra_weight_1, std::mt19937& rng);

    long next(long part) { return buckets[part].next(); }

    long gain_next(long part) { return buckets[part].gain_next(); }

    void move(long v);

    void remove(long v);

    bool done();

    void add_gain(long v, long value);

    bool bucket_done(int part) { return buckets[part].next() == -1; }

  private:
    void init_();

    pmondriaan::hypergraph& H_;
    std::vector<std::vector<long>>& C_;
    std::vector<pmondriaan::gain_buckets> buckets;
    std::vector<long> gains;

    long compute_size_buckets();
};

} // namespace pmondriaan
