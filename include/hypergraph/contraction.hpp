#pragma once

#include <vector>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"

namespace pmondriaan {

/**
 * List of vertices matched.
 */
class match {
  public:
    match(long match, int proc) : id_(match), proc_(proc) {}

    long id() { return id_; }
    int proc() { return proc_; }

  private:
    long id_;
    int proc_;
};

/**
 * A contraction contains the information of how the hypergraph was contracted.
 */
class contraction {
  public:
    contraction() { global_free_weight_ = 0; }

    void add_sample(long id_sample) {
        ids_samples_.push_back(id_sample);
        matches_.push_back(std::vector<pmondriaan::match>());
    }

    void add_samples(pmondriaan::hypergraph& H, const std::vector<long>& indices_samples) {
        for (auto index : indices_samples) {
            ids_samples_.push_back(H(index).id());
        }
        matches_ =
        std::vector<std::vector<pmondriaan::match>>(indices_samples.size(),
                                                    std::vector<pmondriaan::match>());
    }

    void add_match(long sample, long match, int proc) {
        assert(sample >= 0 && (size_t)sample < matches_.size());
        matches_[sample].push_back(pmondriaan::match(match, proc));
    }

    void merge_free_vertices(bulk::world& world, pmondriaan::hypergraph& H);
    void merge_free_vertices(pmondriaan::hypergraph& H);

    /**
     * Assign the free vertices greedily to optimize the weight balance. Returns
     * the weights of the parts.
     */
    std::vector<long> assign_free_vertices(pmondriaan::hypergraph& H,
                                           long max_weight_0,
                                           long max_weight_1,
                                           std::mt19937& rng);

    /**
     * Assign the free vertices greedily to optimize the weight balance. Returns
     * the weights of the parts.
     */
    std::vector<long> assign_free_vertices(bulk::world& world,
                                           pmondriaan::hypergraph& H,
                                           long max_weight_0,
                                           long max_weight_1,
                                           std::mt19937& rng);

    auto& matches(long sample) { return matches_[sample]; }

    long id_sample(long i) { return ids_samples_[i]; }

    auto& free_vertices() { return free_vertices_; }

    long global_free_weight() { return global_free_weight_; }

    size_t size() { return ids_samples_.size(); }

  private:
    std::vector<long> ids_samples_;
    std::vector<std::vector<pmondriaan::match>> matches_;
    std::vector<std::pair<long, long>> free_vertices_;
    long local_free_weight_;
    long global_free_weight_;

    void add_free_vertex_(long id, long weight) {
        free_vertices_.push_back(std::make_pair(id, weight));
    }
    long remove_free_vertices_(pmondriaan::hypergraph& H);
    void assign_all_vertices_(pmondriaan::hypergraph& H, long part);

    /**
     * Assign  the free vertices greedily to optimize the weight balance.
     */
    void assign_free_vertices_(pmondriaan::hypergraph& H,
                               std::vector<long>& weight_parts,
                               long max_weight_0,
                               long max_weight_1,
                               std::mt19937& rng);
};


} // namespace pmondriaan
