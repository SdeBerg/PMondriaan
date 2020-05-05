cmake_minimum_required(VERSION 3.8)
include(FeatureSummary)

project(PMondriaan)

add_subdirectory("ext/bulk")
add_subdirectory("ext/CLI11")
add_subdirectory("ext/googletest")

set(LIB_SOURCES
  "src/hypergraph/readhypergraph.cpp"
  "src/bisect.cpp"
  "src/hypergraph/hypergraph.cpp"
  "src/hypergraph/contraction.cpp"
  "src/recursive_bisection.cpp"
  "src/multilevel_bisect/sample.cpp"
  "src/multilevel_bisect/coarsen.cpp"
  "src/multilevel_bisect/label_propagation.cpp"
  "src/multilevel_bisect/initial_partitioning.cpp"
  "src/multilevel_bisect/uncoarsen.cpp"
  "src/multilevel_bisect/KLFM/KLFM.cpp"
  "src/multilevel_bisect/KLFM/KLFM_parallel.cpp"
  "src/multilevel_bisect/KLFM/gain_buckets.cpp"
  "src/util/write_partitioning.cpp"
)
add_library(PMondriaan ${LIB_SOURCES})

set(EXTERNAL_LIBS
  "bulk_thread"
  "CLI11::CLI11"
)

target_link_libraries(PMondriaan PUBLIC ${EXTERNAL_LIBS})
target_include_directories(PMondriaan PUBLIC "include")

target_compile_features(PMondriaan PUBLIC cxx_std_17)
target_compile_options(PMondriaan PUBLIC
    "-Wfatal-errors"
    "-Werror"
    "-Wall"
    "-O3")

set(TEST_DIR "test")
set(
	TEST_SOURCES
	#"test_bulk.cpp"
	#"test_hypergraph.cpp"
	#"parallel_hypergraph.cpp"
	#"partition.cpp"
	#"random_bisection.cpp"
	#"imbalance.cpp"
	#"redistribute.cpp"
	#"random_sampling.cpp"
	#"renumbering.cpp"
	#"test_CLI11.cpp"
	#"test_label_propagation.cpp"
	#"lp_sampling.cpp"
	#"coarsening.cpp"
	#"bipartition.cpp"
	#"lp_bisect.cpp"
)

foreach(source_file ${TEST_SOURCES})
	string(REPLACE ".cpp" "" source_name ${source_file})
	add_executable(${source_name} ${TEST_DIR}/${source_file})
	target_link_libraries(${source_name} PMondriaan)
endforeach(source_file)

add_executable(Run_PMondriaan "tools/PMondriaan.cpp")
target_link_libraries(Run_PMondriaan PMondriaan)

set(
	UNIT_TEST_SOURCES
	"unittest/test_main.cpp"
	"unittest/hypergraph/readhypergraph_test.cpp"
	"unittest/hypergraph/hypergraph_test.cpp"
	"unittest/hypergraph/contraction_test.cpp"
	"unittest/multilevel_bisect/KLFM/gain_buckets_test.cpp"
	"unittest/multilevel_bisect/KLFM/KLFM_par_test.cpp"
	"unittest/multilevel_bisect/initial_partitioning_test.cpp"
	"unittest/multilevel_bisect/label_propagation_bisect_test.cpp"
	"unittest/multilevel_bisect/bisect_test.cpp"
	"unittest/util/partitioning_to_file_test.cpp"
)

add_executable(pmondriaan_test ${UNIT_TEST_SOURCES})
target_link_libraries(pmondriaan_test PMondriaan gtest)
