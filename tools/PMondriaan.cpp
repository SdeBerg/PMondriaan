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

int main(int argc, char **argv) {
	    
		/* Sequential reading of parameters */
		CLI::App app("PMondriaan settings");
		app.option_defaults()->required();
		
		int p, k;
		double eps, eta;
		std::string matrix_file, hypergraph_weigths, bisection_mode;
		
		app.set_config("--config", "../tools/defaults.toml", "Read a TOML file", true);
		
		auto options = pmondriaan::options();
		CLI::Option *fopt = app.add_option("-f, --file", matrix_file, "File including the hypergraph to be partitioned in matrixmarket format");
		fopt->check(CLI::ExistingFile);
		app.add_option("-k, --k", k, "Number of parts to partition H into");
		app.add_option("-p, --p", p, "Number of processors to be used");
		
		app.add_option("--eps, --epsilon", eps, "Maximum imbalance of the final partitioning");
		app.add_option("--eta", eta, "Maximum imbalance during the parallel computation");
		
		CLI::Option *wopt = app.add_option("--weights", hypergraph_weigths, "How the weights of the vertices should be computed");
		wopt->check(CLI::IsMember({"ones", "degree"}));
		CLI::Option *bopt = app.add_option("--bisect", bisection_mode, "The bisection mode used");
		bopt->check(CLI::IsMember({"random", "multilevel"}));
		CLI::Option *mopt = app.add_option("--metric", options.metric, "Metric to optimized");
		mopt->check(CLI::IsMember({"cutnet", "lambda1"}));
		
		app.add_option("--sample_size", options.sample_size, "The sample size used in the coarsening");
		app.add_option("--max_cluster_size", options.coarsening_max_clustersize, "The maximum clustersize during coarsening");

		CLI11_PARSE(app, argc, argv);
		
		bulk::thread::environment env;
		/* Start parallel part */
		env.spawn(p, [](bulk::world& world) {
			
			
		});

    return 0;
}
