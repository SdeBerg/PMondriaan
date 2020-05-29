#pragma once

namespace pmondriaan {

enum class m : int { cut_net, lambda_minus_one };
enum class bisection : int { random, multilevel };
enum class sampling : int { random, label_propagation };
/**
 *
 */
struct options {
    size_t sample_size;
    size_t coarsening_max_clustersize;
    size_t lp_max_iterations;
    size_t coarsening_nrvertices;
    size_t coarsening_maxrounds;
    size_t KLFM_max_passes;
    size_t KLFM_max_no_gain_moves;
    size_t KLFM_par_number_send_moves;

    m metric;
    bisection bisection_mode;
    sampling sampling_mode;
};

} // namespace pmondriaan
