cmake_minimum_required(VERSION 3.8)
include(FeatureSummary)

project(PMondriaan)

add_subdirectory("ext/bulk/Bulk-develop")
# target_include_directories(bulk INTERFACE "ext/bulk/Bulk-develop")

set(LIB_SOURCES
  "src/hypergraph/readhypergraph.cpp" 
  "src/bisect.cpp"
  "src/hypergraph/hypergraph.cpp"
  "src/recursive_bisection.cpp"
  "src/multilevel_bisect/sample.cpp"
)
add_library(PMondriaan ${LIB_SOURCES})

target_link_libraries(PMondriaan PUBLIC bulk_thread)

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
	"random_sampling.cpp"
)

foreach(source_file ${TEST_SOURCES})
	string(REPLACE ".cpp" "" source_name ${source_file})
	add_executable(${source_name} ${TEST_DIR}/${source_file})
	target_link_libraries(${source_name} PMondriaan)
endforeach(source_file)

