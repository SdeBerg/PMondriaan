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
    match(int match, int proc) : id_(match), proc_(proc) {}

    int id() { return id_; }
    int proc() { return proc_; }

  private:
    int id_;
    int proc_;
};

/**
 * A contraction contains the information of how the hypergraph was contracted.
 */
class contraction {
  public:
    contraction() {
        ids_samples_ = std::vector<int>();
        matches_ = std::vector<std::vector<pmondriaan::match>>();
    }

    void add_samples(pmondriaan::hypergraph& H, std::vector<int> indices_samples) {
        for (auto index : indices_samples) {
            ids_samples_.push_back(H(index).id());
        }
        matches_ =
        std::vector<std::vector<pmondriaan::match>>(indices_samples.size(),
                                                    std::vector<pmondriaan::match>());
    }

    void add_match(int sample, int match, int proc) {
        matches_[sample].push_back(pmondriaan::match(match, proc));
    }

    auto& matches(int sample) { return matches_[sample]; }

    int id_sample(int i) { return ids_samples_[i]; }

    size_t size() { return ids_samples_.size(); }


  private:
    std::vector<int> ids_samples_;
    std::vector<std::vector<pmondriaan::match>> matches_;
};


} // namespace pmondriaan