#include <iostream>
#include <string>

#include <CLI/CLI.hpp>

#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
using environment = bulk::mpi::environment;
#else
#include <bulk/backends/thread/thread.hpp>
using environment = bulk::thread::environment;
#endif

#include <pmondriaan.hpp>

int main(int argc, char** argv) {

    /* Sequential reading of parameters */
    CLI::App app("PMondriaan settings");
    app.option_defaults()->required();

    app.set_config("--config", "../tools/defaults.toml", "Read a TOML file", true);

    struct cli_settings {
        int p = 2;
        int k = 2;
        double eps = 0.05;
        double eta = 0.10;
        std::string matrix_file;
        std::string hypergraph_weights;
        std::string bisection_mode;
        std::string sampling_mode;
        std::string metric;
    };

    cli_settings settings;

    auto options = pmondriaan::options();

    CLI::Option* fopt = app.add_option("-f, --file", settings.matrix_file,
                                       "File including the hypergraph to be "
                                       "partitioned in matrixmarket format");
    fopt->check(CLI::ExistingFile);

    app.add_option("-p, --processors", settings.p,
                   "Number of processors to be used", settings.p);
    app.add_option("-k, --parts", settings.k,
                   "Number of parts to partition H into", settings.k);

    app.add_option("--eps, --epsilon", settings.eps,
                   "Maximum imbalance of the final partitioning", settings.eps);
    app.add_option("--eta", settings.eta,
                   "Maximum imbalance during the parallel computation", settings.eta);

    CLI::Option* wopt =
    app.add_option("--weights", settings.hypergraph_weights,
                   "How the weights of the vertices should be computed");
    CLI::Option* bopt =
    app.add_option("--bisect", settings.bisection_mode, "The bisection mode used");
    CLI::Option* copt =
    app.add_option("--sampling", settings.sampling_mode, "Sampling mode to be used");
    CLI::Option* mopt = app.add_option("--metric", settings.metric, "Metric to optimized");

    wopt->check(CLI::IsMember({"one", "degree"}));
    bopt->check(CLI::IsMember({"random", "multilevel"}));
    copt->check(CLI::IsMember({"random", "label propagation"}));
    mopt->check(CLI::IsMember({"cutnet", "lambda1"}));

    app.add_option("--sample_size", options.sample_size, "The sample size used in the coarsening");
    app.add_option("--max_cluster_size", options.coarsening_max_clustersize,
                   "The maximum clustersize during coarsening");
    app.add_option("--lp_max_iter", options.lp_max_iterations,
                   "The maximum number of iterations in the label propagation");
    app.add_option("--coarsening_nrvertices", options.coarsening_nrvertices,
                   "The number of vertices in the hypergraph when coarsening "
                   "stops");
    app.add_option("--coarsening_max_rounds", options.coarsening_maxrounds,
                   "The maximum number of coarsening rounds");

    CLI11_PARSE(app, argc, argv);

    bulk::thread::environment env;
    /* Start parallel part */
    env.spawn(settings.p, [&settings, &options, &app](bulk::world& world) {
        auto s = world.rank();

        if (s == 0) {
            /* Write the settings to the defaults file */
            auto conf = app.config_to_str();
            std::ofstream out("../tools/settings_run.toml");
            out << conf;
            out.close();
        }

        auto H = pmondriaan::read_hypergraph(settings.matrix_file, world,
                                             settings.hypergraph_weights);

        recursive_bisect(world, H, settings.bisection_mode, settings.sampling_mode,
                         settings.metric, settings.k, settings.eps, settings.eta, options);

        auto lb = pmondriaan::load_balance(world, H, settings.k);
        auto cutsize = pmondriaan::cutsize(world, H, settings.k, settings.metric);


        if (s == 0) {
            world.log("Partitioned hypergraph with %d vertices", H.global_size());
            world.log("Load balance of partitioning found: %lf", lb);
            world.log("Cutsize of partitioning found: %d", cutsize);
        }

        world.sync();
    });

    return 0;
}
