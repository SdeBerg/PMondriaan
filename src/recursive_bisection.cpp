#include <stdlib.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include <iostream>

#include "recursive_bisection.hpp"
#include "hypergraph/hypergraph.hpp"
#include "bisect.hpp"
#include "algorithm.hpp"

namespace pmondriaan {

/**
 * Recursively bisects a hypergraph into k parts with imbalance parameter epsilon.
 * Eta is the load imbalance parameter for the imbalance over the processors during computation.
 */
void recursive_bisect(bulk::world& world, pmondriaan::hypergraph& H, std::string mode, int k, double epsilon, double eta) {
	
	int s = world.rank();
	int p = world.active_processors();
	
	bulk::var<long> local_weight(world);
	local_weight = H.total_weight();
	auto global_weight = bulk::foldl(local_weight, [](auto& lhs, auto rhs) { lhs += rhs; });
	long maxweight = ((1 + epsilon) * global_weight) / k;
	
	// the processors working with processor s on the same part are procs_mypart[0],..., procs_mypart[1] - 1.
	auto procs_mypart = std::vector<int>(2);
	procs_mypart[0] = 0;
	procs_mypart[1] = p;
	
	int start = 0;
	int end = H.size();
	// we need to make sure every processor uses a different seed to generate random numbers
	srand(s + 1);
	int label_low = 0;
	int label_high = k - 1;
	
	// part 0 will always have the heaviest weight
	bulk::var<long> weight_part_0(world);
	bulk::var<long> weight_part_1(world);
	
	int k_, k_low, k_high, q_low, q_high, p_low, p_high, my_part;
	long weight_mypart = global_weight;
	// while we need to give more than one label, we bisect the hypergraph
	while (label_high - label_low > 0) {
		
		k_ = label_high - label_low + 1;
		p = procs_mypart[1] - procs_mypart[0];
		
		double eps =  (double)(maxweight*k_)/(double)weight_mypart - 1.0;
		
		k_low = k_/2;
        k_high = (k_%2==0 ? k_low : k_low+1);
		
		q_low = std::log2(k_low) + 1;
		q_high = std::log2(k_high) + 1;
		
        double delta_low = eps/q_low;
		double delta_high = eps/q_high;
		
        long global_weight_low = ((double)weight_mypart / (double)k_) * (double)k_low * (1.0 + delta_low);
        long global_weight_high = ((double)weight_mypart / (double)k_) * (double)k_high * (1.0 + delta_high);
		
		if (mode == "random") {
			auto weight_parts = bisect_random(world, H, global_weight_high/p, global_weight_low/p, start, end, label_low, label_high);
			weight_part_0 = weight_parts[0];
			weight_part_1 = weight_parts[1];
		}
		
		world.sync();
		
		auto total_weight_0 = pmondriaan::foldl(weight_part_0, [](auto& lhs, auto rhs) { lhs += rhs; }, procs_mypart);
		auto total_weight_1 = pmondriaan::foldl(weight_part_1, [](auto& lhs, auto rhs) { lhs += rhs; }, procs_mypart);
		
		// number of processors working on the low and high part respectively
		p_low = (double)total_weight_1/(double)(total_weight_0 + total_weight_1) * (double)p + 0.5;
		if((k_low == 1) & (k_high > 1)) {
			p_low = 0;
		}
		p_high = p - p_low;
		
		if (s - procs_mypart[0] < p_high){
			my_part = 0;
		}
		else {
			my_part = 1;
		}
		
		//personal low and high label for the next round
		int new_label_low = label_low + my_part * k_high;
		int new_label_high = label_high - (1 - my_part) * k_low;
		
		if (new_label_high - new_label_low > 0) {
			
			long new_max_local_weight; 
			if (my_part == 0) { new_max_local_weight = ceil((double)total_weight_0/(double)p_high) * (1.0 + eta); }
			else { new_max_local_weight = ceil((double)total_weight_1/(double)p_low) * (1.0 + eta); }
			
			/**
			* Redistribute the hypergraph over the processors such that all vertices with label_low are on
			* processors procs_mypart[0] ... procs_mypart[0] + p_high -1 and with label_high on
			* procs_mypart[0] + p_high ... procs_mypart[1] - 1.
			*/
			end = redistribute_hypergraph(world, H, procs_mypart, my_part, label_low, label_high, new_max_local_weight, weight_part_0.value(), weight_part_1.value(), p_low);
			world.log("end %d", end);
			procs_mypart[0] += my_part * p_high;
			procs_mypart[1] -= (1 - my_part) * p_low;
			
			if (my_part == 0) { weight_mypart = total_weight_0; }
			else { weight_mypart = total_weight_1; }
		}
		world.log("weight part %d: %d, weight part %d: %d", label_low, total_weight_0, label_high, total_weight_1);
		
		//personal low and high label for the next round
		label_low = new_label_low;
		label_high = new_label_high;
	}
}

/**
 * Redistributes the hypergraph such that processors with my_part 0 contain all vertices with label_low and
 * all with my_part 1 contain all vertices with label_high. Returns the end of part 0 if part 1 is not assigned
 * any processors.
 */
int redistribute_hypergraph(bulk::world& world, pmondriaan::hypergraph& H, std::vector<int> procs_mypart, int my_part, int label_low, 
					int label_high, long max_local_weight, long weight_part_0, long weight_part_1, int p_low) {
	
	long surplus_0 = (1 - my_part) * max_local_weight - weight_part_0;
	
	auto all_surplus_0 = pmondriaan::gather_all(world, procs_mypart, surplus_0);

	//queue for all received vertices
	auto q = bulk::queue<int, long, int, int[]>(world);
	reduce_surplus(world, H, procs_mypart, label_low, all_surplus_0, q);
	
	//we move all vertices to of part 1 to a separate vector, be be put at the back of the vertices later
	//if part 1 is finished
	auto vertices_1 = std::vector<pmondriaan::vertex>();
	if(p_low > 0) { 
		long surplus_1 = my_part * max_local_weight - weight_part_1;
		auto all_surplus_1 = pmondriaan::gather_all(world, procs_mypart, surplus_1);
		reduce_surplus(world, H, procs_mypart, label_high, all_surplus_1, q);
	}
	else {
		auto i = H.size();
		while(i >= 0) {
			if (H(i).part() == label_low) {
				std::move(H.vertices().begin() + i, H.vertices().begin() + i + 1, std::back_inserter(vertices_1));
				H.vertices().erase(H.vertices().begin() + i);
			}
			else {i--;}
		}
	}
	
	world.sync();
	
	for (auto& [index, weight, part, nets] : q) {
		H.vertices().push_back({index, nets, weight});
		H.vertices().back().set_part(part);
	}	
	
	if(p_low == 0) {
		H.vertices().insert(end(H.vertices()), begin(vertices_1), end(vertices_1));
	}
	return H.size() - vertices_1.size();
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
				q(t).send(H(index).id(), H(index).weight(), H(index).part(), H(index).nets());
				H.vertices().erase(H.vertices().begin() + index);
			}
			
			surplus[s_loc] += total_sent;
			surplus_others = 0;
		}
		t++;
	}
	return 0;
}

} // namespace pmondriaan
