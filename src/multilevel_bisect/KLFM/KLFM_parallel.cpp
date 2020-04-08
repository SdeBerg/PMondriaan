#include <array>
#include <limits>
#include <random>
#include <stack>

#include <bulk/bulk.hpp>
#ifdef BACKEND_MPI
#include <bulk/backends/mpi/mpi.hpp>
#else
#include <bulk/backends/thread/thread.hpp>
#endif

#include "hypergraph/hypergraph.hpp"
#include "multilevel_bisect/KLFM/KLFM_parallel.hpp"
#include "multilevel_bisect/KLFM/gain_buckets.hpp"

namespace pmondriaan {

/**
 * Runs the KLFM algorithm in parallel to improve a given partitioning. Return the cutsize of the best solution.
 */
long KLFM_par(bulk::world& world,
              pmondriaan::hypergraph& H,
              std::vector<std::vector<long>>& C,
              long weight_0,
              long weight_1,
              long max_weight_0,
              long max_weight_1,
              pmondriaan::options& opts,
              std::mt19937& rng,
              long cut_size) {

    int pass = 0;
    long prev_cut_size;

    // TODO: use C to compute cutsize without communication
    if (cut_size == std::numeric_limits<decltype(cut_size)>::max()) {
        prev_cut_size = pmondriaan::cutsize(world, H, opts.metric);
    } else {
        prev_cut_size = cut_size;
    }
    auto total_weights = std::array<long, 2>();
    total_weights[0] = weight_0;
    total_weights[1] = weight_1;

    while (pass < opts.KLFM_max_passes) {
        auto result = KLFM_pass_par(world, H, C, prev_cut_size, total_weights,
                                    max_weight_0, max_weight_1, opts, rng);
        if (result < prev_cut_size) {
            prev_cut_size = result;
        } else {
            break;
        }
        pass++;
    }
    return prev_cut_size;
}

/**
 * Runs a single pass of the KLFM algorithm in parallel to improve a given partitioning.
 */
long KLFM_pass_par(bulk::world& world,
                   pmondriaan::hypergraph& H,
                   std::vector<std::vector<long>>& C,
                   long cut_size,
                   std::array<long, 2>& total_weights,
                   long max_weight_0,
                   long max_weight_1,
                   pmondriaan::options& opts,
                   std::mt19937& rng) {
    int s = world.rank();
    int p = world.active_processors();

    auto gain_structure = pmondriaan::gain_structure(H, C);

    long best_cut_size = cut_size;
    auto no_improvement_moves = std::vector<int>();

    // Each processor is responsible for keeping track of some nets
    auto net_partition = bulk::block_partitioning<1>({H.global_number_nets()}, {p});
    auto previous_C = std::vector<std::vector<long>>(net_partition.local_count(s),
                                                     std::vector<long>(2));
    for (auto i = 0; i < net_partition.local_count(s); i++) {
        previous_C[i][0] = C[net_partition.global(i, s)[0]][0];
        previous_C[i][1] = C[net_partition.global(i, s)[0]][1];
    }

    long cut_size_my_nets = 0;
    for (auto i = 0; i < net_partition.local_count(s); i++) {
        auto n = net_partition.global(i, s)[0];
        if ((C[n][0] == 0) || (C[n][1] == 0)) {
            cut_size_my_nets += H.net(n).cost();
        }
    }

    auto moves_queue = bulk::queue<long, long, int>(world);
    auto update_nets = bulk::queue<int, long>(world);
    bulk::var<int> rejected(world);
    bulk::var<long> new_weight_0(world);
    bulk::var<long> new_weight_1(world);
    bulk::var<bool> done(world);

    bool all_done = false;
    while (!all_done) {
        auto prev_total_weights = total_weights;
        // Find best KLFM_par_number_send_moves moves
        auto moves =
        std::vector<std::tuple<int, long, long>>(opts.KLFM_par_number_send_moves);
        find_top_moves(H, gain_structure, moves, total_weights, max_weight_0,
                       max_weight_1, rng);
        world.sync();

        // Send moves to processor 0
        for (auto& move : moves) {
            moves_queue(0).send(std::get<1>(move), std::get<2>(move), s);
        }
        world.sync();

        if (s == 0) {
            auto rejected_all = reject_unbalanced_moves(p, moves_queue, prev_total_weights,
                                                        max_weight_0, max_weight_1);
            for (auto t = 0; t < p; t++) {
                rejected(t) = rejected_all[t];
            }
            new_weight_0.broadcast(prev_total_weights[0]);
            new_weight_1.broadcast(prev_total_weights[1]);
        }
        world.sync();
        total_weights[0] = new_weight_0.value();
        total_weights[1] = new_weight_1.value();

        // We reverse moves that have not been selected
        auto index = moves.size() - 1;
        if (rejected > 0) {
            while (rejected != 0) {
                if (std::get<2>(moves[index]) > 0) {
                    auto& vertex = H(H.local_id(std::get<0>(moves[index])));
                    H.move(vertex.id(), C);
                    rejected--;
                }
                index--;
            }
        } else {
            while (rejected != 0) {
                if (std::get<2>(moves[index]) < 0) {
                    auto& vertex = H(H.local_id(std::get<0>(moves[index])));
                    H.move(vertex.id(), C);
                    rejected++;
                }
                index--;
            }
        }

        // We send updates about counts in nets to responsible processors
        // TODO: Make this more efficient by not sending everything
        for (auto i = 0u; i < H.nets().size(); i++) {
            update_nets(net_partition.owner(i)).send(i, C[i][0]);
        }
        world.sync();

        // Update the counts of my nets
        for (auto i = 0; i < net_partition.local_count(s); i++) {
            C[net_partition.global(i, s)[0]][0] = previous_C[i][0];
            C[net_partition.global(i, s)[0]][1] = previous_C[i][1];
        }
        for (const auto& [net, C_0] : update_nets) {
            auto change = C_0 - previous_C[net_partition.local(net)[0]][0];
            if (change != 0) {
                if ((C[net][0] == 0) || (C[net][1] == 0)) {
                    cut_size_my_nets += H.net(net).cost();
                }
                C[net][0] += change;
                C[net][1] -= change;
                // TODO make sure we also adjust the gains here
                if ((C[net][0] == 0) || (C[net][1] == 0)) {
                    cut_size_my_nets -= H.net(net).cost();
                }
            }
        }
        world.sync();
        // Send updated counts to all processors
        // TODO: improve to send only to processors that need it
        for (auto i = 0; i < net_partition.local_count(s); i++) {
            auto global_id = net_partition.global(i, s)[0];
            for (auto t = 0; t < p; t++) {
                update_nets(t).send(i, C[i][0]);
            }
            previous_C[i][0] = C[global_id][0];
            previous_C[i][1] = C[global_id][1];
        }
        world.sync();

        // TODO: Update C and update gains
        for (const auto& [net, C_0] : update_nets) {
            if (C_0 != C[net][0]) {
                if ((C[net][0] == 0) || (C[net][1] == 0)) {
                }
                C[net][0] = C_0;
                C[net][1] = H.net(net).size() - C_0;
            }
        }

        // We also send all processors the total cutsize of the nets this p is responsible for
        cut_size = bulk::sum(world, cut_size_my_nets);

        if (cut_size > best_cut_size) {
            for (auto& move : moves) {
                no_improvement_moves.push_back(std::get<0>(move));
            }
        } else {
            best_cut_size = cut_size;
            no_improvement_moves.clear();
        }

        if (gain_structure.done()) {
            done = true;
        }
        all_done =
        bulk::foldl(done, [](auto& lhs, auto rhs) { lhs = (lhs && rhs); });
        all_done = true;
    }

    // Reverse moves up to best point
    long weight_change = 0;
    for (auto v : no_improvement_moves) {
        auto& vertex = H(H.local_id(v));
        if (vertex.part() == 0) {
            weight_change -= vertex.weight();
        } else {
            weight_change += vertex.weight();
        }
        H.move(v, C);
    }
    auto total_change = bulk::sum(world, weight_change);
    total_weights[0] += total_change;
    total_weights[1] -= total_change;

    return best_cut_size;
}

/**
 * Finds the best moves for a processor sequentially, by only updating local data.
 */
void find_top_moves(pmondriaan::hypergraph& H,
                    pmondriaan::gain_structure& gain_structure,
                    std::vector<std::tuple<int, long, long>>& moves,
                    std::array<long, 2>& weights,
                    long max_weight_0,
                    long max_weight_1,
                    std::mt19937& rng) {
    auto max_extra_weight = std::array<long, 2>();
    auto moves_found = 0u;

    while (!gain_structure.done() && (moves_found < moves.size())) {
        max_extra_weight[0] = max_weight_0 - weights[0];
        max_extra_weight[1] = max_weight_1 - weights[1];

        int part_to_move =
        gain_structure.part_next(max_extra_weight[0], max_extra_weight[1], rng);
        auto v_to_move = gain_structure.next(part_to_move);

        if (max_extra_weight[(part_to_move + 1) % 2] - H(H.local_id(v_to_move)).weight() >= 0) {
            auto weight_v = H(H.local_id(v_to_move)).weight();
            if (part_to_move == 0) {
                moves[moves_found] =
                std::make_tuple(v_to_move, gain_structure.gain_next(part_to_move),
                                -1 * weight_v);
                weights[0] -= weight_v;
                weights[1] += weight_v;
            } else {
                moves[moves_found] =
                std::make_tuple(v_to_move, gain_structure.gain_next(part_to_move), weight_v);
                weights[1] -= weight_v;
                weights[0] += weight_v;
            }
            gain_structure.move(v_to_move);
            moves_found++;
        } else {
            gain_structure.remove(v_to_move);
        }
    }
}

/**
 * Determines how many moves from part 0 or 1 should be rejected on each processor.
 * Return the number of rejected moves for each processor. This is positive when
 * we have to move back vertices from part 0 to part 1 and positive otherwise.
 */
std::vector<int> reject_unbalanced_moves(int p,
                                         bulk::queue<long, long, int>& moves_queue,
                                         std::array<long, 2>& total_weights,
                                         long max_weight_0,
                                         long max_weight_1) {
    // The number of moves to be rejected for each processor
    auto rejected = std::vector<int>(p, 0);
    auto done = std::vector<bool>(p, false);

    // Sort queue on gain
    std::sort(moves_queue.begin(), moves_queue.end());
    std::stack<std::pair<long, int>> received_moves;
    long total_balance = 0;
    for (const auto& [gain, weight_change, t] : moves_queue) {
        total_balance += weight_change;
        received_moves.push(std::make_pair(weight_change, t));
    }
    total_weights[0] += total_balance;
    total_weights[1] -= total_balance;
    while (total_weights[0] > max_weight_0) {
        if (received_moves.empty()) {
            break;
        }
        if ((received_moves.top().first < 0) && (!done[received_moves.top().second])) {
            if (total_weights[1] - received_moves.top().first <= max_weight_1) {
                rejected[received_moves.top().second]++;
                total_weights[0] += received_moves.top().first;
                total_weights[1] -= received_moves.top().first;
            } else {
                done[received_moves.top().second] = true;
            }
        }
        received_moves.pop();
    }

    while (total_weights[1] > max_weight_1) {
        if (received_moves.empty()) {
            break;
        }
        if ((received_moves.top().first > 0) && (!done[received_moves.top().second])) {
            if (total_weights[0] + received_moves.top().first <= max_weight_0) {
                rejected[received_moves.top().second]--;
                total_weights[0] += received_moves.top().first;
                total_weights[1] -= received_moves.top().first;
            } else {
                done[received_moves.top().second] = true;
            }
        }
        received_moves.pop();
    }
    return rejected;


    /*for (const auto& [gain, weight_change, t] : moves_queue) {
        if (stops[t] == -1) {
            auto new_w_0 = total_weights[0] + weight_change;
            auto new_w_1 = total_weights[1] - weight_change;
            if ((new_w_0 <= max_weight_0) && (new_w_1 <= max_weight_1)) {
                current[t]++;
                total_weights[0] = new_w_0;
                total_weights[1] = new_w_1;
            } else {
                number_done++;
                stops[t] = current[t];
                if (number_done == world.active_processors()) {
                    break;
                }
            }
        }
    }

    for (auto t = 0; t < world.active_processors(); t++) {
        if (stops[t] == -1) {
            stop(t) = current[t];
        } else {
            stop(t) = stops[t];
        }
    }*/
}

} // namespace pmondriaan
