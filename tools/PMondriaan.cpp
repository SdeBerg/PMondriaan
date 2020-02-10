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
	
	app.set_config("--config", "../tools/defaults.toml", "Read a TOML file", true);
	
	int p, k;
	double eps, eta;
	std::string matrix_file, hypergraph_weigths, bisection_mode, sampling_mode, metric;
	auto options = pmondriaan::options();
	
	CLI::Option *fopt = app.add_option("-f, --file", matrix_file, "File including the hypergraph to be partitioned in matrixmarket format");
	fopt->check(CLI::ExistingFile);
	
	app.add_option("-p, --p", p, "Number of processors to be used");
	app.add_option("-k, --k", k, "Number of parts to partition H into");
	
	app.add_option("--eps, --epsilon", eps, "Maximum imbalance of the final partitioning");
	app.add_option("--eta", eta, "Maximum imbalance during the parallel computation");
	
	CLI::Option *wopt = app.add_option("--weights", hypergraph_weigths, "How the weights of the vertices should be computed");
	CLI::Option *bopt = app.add_option("--bisect", bisection_mode, "The bisection mode used");
	CLI::Option *copt = app.add_option("--sampling", sampling_mode, "Sampling mode to be used");
	CLI::Option *mopt = app.add_option("--metric", metric, "Metric to optimized");
	
	wopt->check(CLI::IsMember({"one", "degree"}));
	bopt->check(CLI::IsMember({"random", "multilevel"}));
	copt->check(CLI::IsMember({"random", "label propagation"}));
	mopt->check(CLI::IsMember({"cutnet", "lambda1"}));
	
	app.add_option("--sample_size", options.sample_size, "The sample size used in the coarsening");
	app.add_option("--max_cluster_size", options.coarsening_max_clustersize, "The maximum clustersize during coarsening");
	app.add_option("--lp_max_iter", options.lp_max_iterations, "The maximum number of iterations in the label propagation");

	CLI11_PARSE(app, argc, argv);
	
	/* Write the settings to the defaults file */
	auto conf = app.config_to_str();
	std::ofstream out("../tools/defaults.toml");
	out << conf;
	out.close();
	
	bulk::thread::environment env;
	/* Start parallel part */
	env.spawn(p, [](bulk::world& world) {
		auto p = world.active_processors();
		auto s = world.rank();
		
		auto q_settings = bulk::queue<int, double, double, std::string, std::string, std::string, std::string, std::string, pmondriaan::options>(world);

		int k;
		double eps, eta;
		std::string matrix_file, hypergraph_weigths, bisection_mode, sampling_mode, metric;

		auto options = pmondriaan::options();
		
		if (s == 0) {
			
			/* Sequential reading of parameters */
			CLI::App app("PMondriaan settings");
			app.prefix_command(true);
			app.option_defaults()->required();
			
			app.set_config("--config", "../tools/defaults.toml", "Read a TOML file", true);
			
			app.add_option("-f, --file", matrix_file, "File including the hypergraph to be partitioned in matrixmarket format");
			
			app.add_option("-k, --k", k, "Number of parts to partition H into");
			app.add_option("--eps, --epsilon", eps, "Maximum imbalance of the final partitioning");
			app.add_option("--eta", eta, "Maximum imbalance during the parallel computation");
			app.add_option("--weights", hypergraph_weigths, "How the weights of the vertices should be computed");
			app.add_option("--bisect", bisection_mode, "The bisection mode used");
			app.add_option("--sampling", sampling_mode, "Sampling mode to be used");
			app.add_option("--metric", metric, "Metric to optimized");
			app.add_option("--sample_size", options.sample_size, "The sample size used in the coarsening");
			app.add_option("--max_cluster_size", options.coarsening_max_clustersize, "The maximum clustersize during coarsening");
			app.add_option("--lp_max_iter", options.lp_max_iterations, "The maximum number of iterations in the label propagation");

			const char  arg0[] = "empty";
			const char* argv[] = { &arg0[0], NULL };
			int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
			app.parse(argc, &argv[0]);
			
			for (int t = 0; t < p; t++) {
				q_settings(t).send(k, eps, eta, matrix_file, hypergraph_weigths, bisection_mode, sampling_mode, metric, options);
			}
			
		}
		
		world.sync();
		
		for (const auto& [k_, eps_, eta_, matrix_file_, hypergraph_weigths_, bisection_mode_, sampling_mode_, metric_, options_] : q_settings) {
            k = k_;
			eps = eps_;
			eta = eta_;
			matrix_file = matrix_file_;
			hypergraph_weigths = hypergraph_weigths_;
			bisection_mode = bisection_mode_;
			sampling_mode = sampling_mode_;
			metric = metric_;
			options = options_;
        }
		
		auto H = pmondriaan::read_hypergraph(matrix_file, world, hypergraph_weigths);

		recursive_bisect(world, H, bisection_mode, sampling_mode, metric, k, eps, eta, options);
		
	});

    return 0;
}
