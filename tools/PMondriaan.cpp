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
		std::string matrix_file, hypergraph_weigths, bisection_mode;
		auto options = pmondriaan::options();
		
		CLI::Option *fopt = app.add_option("-f, --file", matrix_file, "File including the hypergraph to be partitioned in matrixmarket format");
		fopt->check(CLI::ExistingFile);
		
		app.add_option("-p, --p", p, "Number of processors to be used");
		app.add_option("-k, --k", k, "Number of parts to partition H into");
		
		app.add_option("--eps, --epsilon", eps, "Maximum imbalance of the final partitioning");
		app.add_option("--eta", eta, "Maximum imbalance during the parallel computation");
		
		CLI::Option *wopt = app.add_option("--weights", hypergraph_weigths, "How the weights of the vertices should be computed");
		CLI::Option *bopt = app.add_option("--bisect", bisection_mode, "The bisection mode used");
		CLI::Option *mopt = app.add_option("--metric", options.metric, "Metric to optimized");
		
		wopt->check(CLI::IsMember({"ones", "degree"}));
		bopt->check(CLI::IsMember({"random", "multilevel"}));
		mopt->check(CLI::IsMember({"cutnet", "lambda1"}));
		
		app.add_option("--sample_size", options.sample_size, "The sample size used in the coarsening");
		app.add_option("--max_cluster_size", options.coarsening_max_clustersize, "The maximum clustersize during coarsening");

		CLI11_PARSE(app, argc, argv);
		
		auto conf = app.config_to_str();
		std::ofstream out("../tools/defaults.toml");
		out << conf;
		out.close();
		
		bulk::thread::environment env;
		/* Start parallel part */
		env.spawn(p, [](bulk::world& world) {
			//auto p = world.active_processors();
			auto s = world.rank();
			
			if (s == 0) {
				
				/* Sequential reading of parameters */
				CLI::App app("PMondriaan settings");
				app.prefix_command(true);
				app.option_defaults()->required();
				
				app.set_config("--config", "../tools/defaults.toml", "Read a TOML file", true);
				
				int k;
				double eps, eta;
				std::string matrix_file, hypergraph_weigths, bisection_mode;

				auto options = pmondriaan::options();
				
				app.add_option("-f, --file", matrix_file, "File including the hypergraph to be partitioned in matrixmarket format");
				
				app.add_option("-k, --k", k, "Number of parts to partition H into");
				app.add_option("--eps, --epsilon", eps, "Maximum imbalance of the final partitioning");
				app.add_option("--eta", eta, "Maximum imbalance during the parallel computation");
				app.add_option("--weights", hypergraph_weigths, "How the weights of the vertices should be computed");
				app.add_option("--bisect", bisection_mode, "The bisection mode used");
				app.add_option("--metric", options.metric, "Metric to optimized");
				app.add_option("--sample_size", options.sample_size, "The sample size used in the coarsening");
				app.add_option("--max_cluster_size", options.coarsening_max_clustersize, "The maximum clustersize during coarsening");

				const char  arg0[] = "empty";
				const char* argv[] = { &arg0[0], NULL };
				int argc = (int)(sizeof(argv) / sizeof(argv[0])) - 1;
				app.parse(argc, &argv[0]);
				
				world.log("k: %d", k);
				
			}
			/*
			
			auto H = pmondriaan::read_hypergraph(matrix_file, world, hypergraph_weigths);

			recursive_bisect(world, H, bisection_mode, k, eps, eta);*/
			
		});

    return 0;
}
