#pragma once

namespace pmondriaan {

enum class m : long { cut_net, lambda_minus_one };
enum class bisection : long { random, multilevel };
enum class sampling : long { random, label_propagation };
/**
 *
 */
class options {
  public:
    long sample_size;
    long coarsening_max_clustersize;
    long lp_max_iterations;
    long coarsening_nrvertices;
    long coarsening_maxrounds;
    long KLFM_max_passes;
    long KLFM_par_number_send_moves;

    m metric;
    bisection bisection_mode;
    sampling sampling_mode;
};

} // namespace pmondriaan