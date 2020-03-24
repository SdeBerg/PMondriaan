#pragma once

namespace pmondriaan {

enum class m : int { cut_net, lambda_minus_one };
enum class bisection : int { random, multilevel };
enum class sampling : int { random, label_propagation };
/**
 *
 */
class options {
  public:
    int sample_size;
    int coarsening_max_clustersize;
    int lp_max_iterations;
    int coarsening_nrvertices;
    int coarsening_maxrounds;
    int KLFM_max_passes;

    m metric;
    bisection bisection_mode;
    sampling sampling_mode;
};

} // namespace pmondriaan