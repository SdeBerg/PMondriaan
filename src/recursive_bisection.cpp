#include <stdlib.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <stack>

#include <iostream>

#include "recursive_bisection.hpp"
#include "hypergraph/hypergraph.hpp"
#include "bisect.hpp"
#include "algorithm.hpp"
#include "work_item.hpp"
#include "options.hpp"

namespace pmondriaan {

/**
 * Recursively bisects a hypergraph into k parts with imbalance parameter epsilon.
 * Eta is the load imbalance parameter for the imbalance over the processors during computation.
 */
void recursive_bisect(bulk::world& world, pmondriaan::hypergraph& H, std::string mode, std::string sampling_mode, std::string metric, 
				int k, double epsilon, double eta, pmondriaan::options& opts) {
	
	int s = world.rank();
	int p = world.active_processors();
	
	// we need to make sure every processor uses a different seed to generate random numbers
	srand(s + 1);
	
	//we first compute the max weight of each part
	bulk::var<long> local_weight(world);
	local_weight = H.total_weight();
	auto global_weight = bulk::foldl(local_weight, [](auto& lhs, auto rhs) { lhs += rhs; });
	long maxweight = ((1 + epsilon) * global_weight) / k;
	
	// the processors working with processor s on the same part are procs_mypart[0],..., procs_mypart[1] - 1.
	auto procs_mypart = std::vector<int>(2);
	procs_mypart[0] = 0;
	procs_mypart[1] = p;
	
	//start and end index of the current part to be split
	int start = 0;
	int end = H.size();
	auto jobs = std::stack<pmondriaan::work_item>();

	int label_low = 0;
	int label_high = k - 1;
	
	// part 0 will always have the smallest weight
	bulk::var<long> weight_part_0(world);
	bulk::var<long> weight_part_1(world);
	
	int k_, k_low, k_high, p_low, p_high, my_part;
	long weight_mypart = global_weight;
	
	if (p == 1) {
		jobs.push(pmondriaan::work_item(start, end, label_low, label_high, weight_mypart));
	}
	
	// while we need to give more than one label and need to use more than one processor, we bisect the hypergraph in parallel
	while (label_high - label_low > 0 && p > 1) {
		
		k_ = label_high - label_low + 1;
		p = procs_mypart[1] - procs_mypart[0];
		
		k_low = k_/2;
        k_high = (k_%2==0 ? k_low : k_low+1);
		
		auto max_global_weights = compute_max_global_weight(k_, k_low, k_high, weight_mypart, maxweight);
		
		if (mode == "random") {
			auto weight_parts = bisect_random(world, H, max_global_weights[0]/p, max_global_weights[1]/p, start, end, label_low, label_high);
			weight_part_0 = weight_parts[0];
			weight_part_1 = weight_parts[1];
		}
		
		if (mode == "multilevel") {
			auto weight_parts = bisect_multilevel(world, H, opts, sampling_mode, metric, max_global_weights[0], max_global_weights[1], start, end, label_low, label_high);
			weight_part_0 = weight_parts[0];
			weight_part_1 = weight_parts[1];
		}
		
		world.sync();
		
		auto total_weight_0 = pmondriaan::foldl(weight_part_0, [](auto& lhs, auto rhs) { lhs += rhs; }, procs_mypart);
		auto total_weight_1 = pmondriaan::foldl(weight_part_1, [](auto& lhs, auto rhs) { lhs += rhs; }, procs_mypart);

		// number of processors working on the low and high part respectively
		p_low = (double)total_weight_0/(double)(total_weight_0 + total_weight_1) * (double)p + 0.5;
		if((k_low == 1) & (k_high > 1)) {
			p_low = 0;
		}
		p_high = p - p_low;
		
		if (s - procs_mypart[0] < p_low){
			my_part = 0;
		}
		else {
			my_part = 1;
		}
		
		//personal low and high label for the next round
		int new_label_low = label_low + my_part * k_low;
		int new_label_high = label_high - (1 - my_part) * k_high;

		if (new_label_high - new_label_low > 0) {
			
			long new_max_local_weight; 
			if (my_part == 0) { new_max_local_weight = ceil((double)total_weight_0/(double)p_low) * (1.0 + eta); }
			else { new_max_local_weight = ceil((double)total_weight_1/(double)p_high) * (1.0 + eta); }
			
			/**
			* Redistribute the hypergraph over the processors such that all vertices with label_low are on
			* processors procs_mypart[0] ... procs_mypart[0] + p_high -1 and with label_high on
			* procs_mypart[0] + p_high ... procs_mypart[1] - 1.
			*/		
			end = redistribute_hypergraph(world, H, procs_mypart, my_part, label_low, label_high, new_max_local_weight, weight_part_0.value(), weight_part_1.value(), p_low);

			procs_mypart[0] += my_part * p_low;
			procs_mypart[1] -= (1 - my_part) * p_high;

			p = procs_mypart[1] - procs_mypart[0];
			
			if (my_part == 0) { weight_mypart = total_weight_0; }
			else { weight_mypart = total_weight_1; }
			
			//if we are the only processor left, we add a sequential job
			if (p == 1) {
				jobs.push(pmondriaan::work_item(start, end, new_label_low, new_label_high, weight_mypart));
			}
		}
		
		world.log("weight part %d: %d, weight part %d: %d", label_low, total_weight_0, label_high, total_weight_1);
		
		//personal low and high label for the next round
		label_low = new_label_low;
		label_high = new_label_high;
	}
	
	//we do the rest of the work sequentially
	while (!jobs.empty()) {
		auto job = jobs.top();
		jobs.pop();
		
		start = job.start();
		end = job.end();
		label_low = job.label_low();
		label_high = job.label_high();
		weight_mypart = job.weight();
		
		while(label_high - label_low > 0) {
			k_ = label_high - label_low + 1;
			
			k_low = k_/2;
			k_high = (k_%2==0 ? k_low : k_low+1);
			
			auto max_global_weights = compute_max_global_weight(k_, k_low, k_high, weight_mypart, maxweight);
			
			if (mode == "random") {
				auto weight_parts = bisect_random(world, H, max_global_weights[0], max_global_weights[1], start, end, label_low, label_high);
				weight_part_0 = weight_parts[0];
				weight_part_1 = weight_parts[1];
			}
			
			if (mode == "multilevel") {
				auto weight_parts = bisect_multilevel(world, H, opts, sampling_mode, metric, max_global_weights[0], max_global_weights[1], start, end, label_low, label_high);
				weight_part_0 = weight_parts[0];
				weight_part_1 = weight_parts[1];
			}
			
			int label_low_0 = label_low;
			int label_high_0 = label_high - k_high;
			int label_low_1 = label_low + k_low;
			int label_high_1 = label_high;
			
			int end_1 = end;
			
			if (label_high_0 - label_low_0 > 0) {
				reorder_hypergraph(H, start, end, label_low, label_high);
			}
			if (label_high_1 - label_low_1 > 0) {
				jobs.push(pmondriaan::work_item(end, end_1, label_low_1, label_high_1, weight_part_1.value()));
			}
			world.log("weight part %d: %d, weight part %d: %d", label_low, weight_part_0.value(), label_high, weight_part_1.value());
			label_low = label_low_0;
			label_high = label_high_0;
		}
	}
	
	world.sync();
}

std::vector<long> compute_max_global_weight(int k_, int k_low, int k_high, long weight_mypart, long maxweight) {
	
		double eps =  (double)(maxweight*k_)/(double)weight_mypart - 1.0;
		
		int q_low = std::log2(k_low) + 1;
		int q_high = std::log2(k_high) + 1;
		
        double delta_low = eps/q_low;
		double delta_high = eps/q_high;
		
        long global_weight_low = ((double)weight_mypart / (double)k_) * (double)k_low * (1.0 + delta_low);
        long global_weight_high = ((double)weight_mypart / (double)k_) * (double)k_high * (1.0 + delta_high);

		return {global_weight_low, global_weight_high};
}

/**
 * Redistributes the hypergraph such that processors with my_part 0 contain all vertices with label_low and
 * all with my_part 1 contain all vertices with label_high. Returns the end of part 0 if part 1 is not assigned
 * any processors.
 */
int redistribute_hypergraph(bulk::world& world, pmondriaan::hypergraph& H, std::vector<int> procs_mypart, int my_part, int label_low, 
					int label_high, long max_local_weight, long weight_part_0, long weight_part_1, int p_low) {
	
	long surplus_0 = (1 - my_part) * max_local_weight - weight_part_0;
	long surplus_1 = my_part * max_local_weight - weight_part_1;
	auto all_surplus_0 = pmondriaan::gather_all(world, procs_mypart, surplus_0);
	auto all_surplus_1 = pmondriaan::gather_all(world, procs_mypart, surplus_1);
	
	//we move all vertices to of part 1 to a separate vector, be be put at the back of the vertices later
	//if part 1 is finished
	auto vertices_0 = std::vector<pmondriaan::vertex>();
	
	//queue for all received vertices
	auto q = bulk::queue<int, long, int, int[]>(world);
	reduce_surplus(world, H, procs_mypart, label_high, all_surplus_1, q);
	
	if(p_low > 0) {
		reduce_surplus(world, H, procs_mypart, label_low, all_surplus_0, q);
	}
	else {
		long i = H.size();
		while(i >= 0) {
			if (H(i).part() == label_low) {
				std::move(H.vertices().begin() + i, H.vertices().begin() + i + 1, std::back_inserter(vertices_0));
				H.vertices().erase(H.vertices().begin() + i);
			}
			i--;
		}
	}
	
	world.sync();

	for (auto& [id, weight, part, nets] : q) {
		H.vertices().push_back({id, nets, weight});
		H.vertices().back().set_part(part);
		H.add_to_nets(H.vertices().back());
	}	
	
	
	if(p_low == 0) {
		H.vertices().insert(end(H.vertices()), begin(vertices_0), end(vertices_0));
	}
	
	H.renumber_vertices();
	
	return H.size() - vertices_0.size();
}

/**
 * Distributes the surplus vertices of a processor over the other processors. Returns 0 if
 * all surplus is gone, -1 if there is still surplus left.
 */		
int reduce_surplus(bulk::world& world, pmondriaan::hypergraph& H, std::vector<int> procs_mypart, 
				int label, bulk::coarray<long>& surplus, bulk::queue<int, long, int, int[]>& q) {

	int s_loc = world.rank() - procs_mypart[0];
	int p_loc = procs_mypart[1] - procs_mypart[0];

	if (surplus[s_loc] >= 0) {
		return 0;
	}
	
	//index of next vertex to be sent
	int index = 0;
	long surplus_others = 0;
	int t = s_loc + 1;
	while (surplus[s_loc] < 0) {
		if (t >= p_loc) {
			t = 0;
		}
		if (t == s_loc){
			std::cerr << "Error: failed to lose all surplus\n";
			return -1;
		}
		surplus_others += surplus[t];
		if (surplus[t] > 0 && surplus_others > 0) {

			long max = std::min(surplus[t], surplus_others);
			max = std::min(max, -1*surplus[s_loc]);
			
			long total_sent = 0;
			while (total_sent < max) {
				//search for the next vertex to be sent
				while(H(index).part() != label) {
					index++;
				}
				
				total_sent += H(index).weight();
				q(t + procs_mypart[0]).send(H(index).id(), H(index).weight(), H(index).part(), H(index).nets());
				
				int id = H(index).id();
				H.remove_from_nets(id);
				H.map().erase(id);
				H.vertices().erase(H.vertices().begin() + index);
			}
			
			surplus[s_loc] += total_sent;
			surplus_others = 0;
		}
		t++;
	}
	
	return 0;
}

void reorder_hypergraph(pmondriaan::hypergraph& H, int start, int& end, int label_low, int label_high) {
	int pivot = start;
	while (pivot != end - 1) {
		if (H(pivot).part() == label_high) {
			while ((H(end - 1).part() != label_low) && (end - 1 != pivot)) {
				end--;
			}
			std::swap(H(pivot),H(end - 1));
		}
		pivot++;
	}
}

} // namespace pmondriaan
