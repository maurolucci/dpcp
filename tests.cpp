#include "bp.hpp"
#include "col.hpp"
#include "compact_ilp.hpp"
#include "graph.hpp"
#include "heur.hpp"
#include "lp.hpp"
#include "params.hpp"
#include "random.hpp"
#include "stats.hpp"

#include <algorithm> // shuffle
#include <filesystem>
#include <iostream>
#include <random> // std::default_random_engine

using recursive_directory_iterator =
    std::filesystem::recursive_directory_iterator;

std::string SEPBAR =
    "*************************************************************";

int main() {

  // Set seed
  rng.seed(0);

  for (const auto &file : recursive_directory_iterator("input/tests/")) {

    if (file.path().extension() != ".txt")
      continue;

    std::cout << "*** Solving: " << file << " ***" << std::endl;

    // Open file
    std::ifstream in(file.path());

    // Read hypergraph
    HGraph hg;
    read_hypergrah(hg, in);
    std::cout << "Hypergraph:" << std::endl;
    std::cout << "\tVertices: " << hg.nbVertices() << std::endl;
    std::cout << "\tHyperedges: " << hg.hyperedges().size() << std::endl;

    // Read parameters
    Params params;

    // Build conflict graph
    Graph graph;
    get_conflict_graph(hg, graph);
    std::cout << "Conflict graph:" << std::endl;
    std::cout << "\tVertices: " << num_vertices(graph) << std::endl;
    std::cout << "\tEdges: " << num_edges(graph) << std::endl;

    std::cout << std::endl << SEPBAR << std::endl << std::endl;

    // Solve with one-step heuristic
    GraphEnv genv(graph, params);
    Col col0s;
    std::cout << "1. Running two-step heuristic for DPCP..." << std::endl;
    Stats stats0s = dpcp_1_step_greedy_heur(genv, col0s);
    std::cout << "Time: " << stats0s.time << std::endl;
    if (stats0s.state == FEASIBLE)
      std::cout << "Value: " << stats0s.ub << std::endl;
    else
      std::cout << "Value: No solution found :(" << std::endl;

    std::cout << std::endl << SEPBAR << std::endl << std::endl;

    // Solve with two-step greedy heuristic
    Col col1s;
    std::cout << "2. Running two-step greedy heuristic for DPCP..."
              << std::endl;
    Stats stats1s = dpcp_2_step_greedy_heur(genv, col1s);
    std::cout << "Time: " << stats1s.time << std::endl;
    if (stats1s.state == FEASIBLE)
      std::cout << "Value: " << stats1s.ub << std::endl;
    else
      std::cout << "Value: No solution found :(" << std::endl;

    std::cout << std::endl << SEPBAR << std::endl << std::endl;

    // Solve with two-step semigreedy heuristic
    Col col2s;
    std::cout << "2. Running two-step semigreedy heuristic for DPCP..."
              << std::endl;
    Stats stats2s = dpcp_2_step_semigreedy_heur(genv, col2s, 100);
    std::cout << "Time: " << stats2s.time << std::endl;
    if (stats2s.state == FEASIBLE)
      std::cout << "Value: " << stats2s.ub << std::endl;
    else
      std::cout << "Value: No solution found :(" << std::endl;

    std::cout << std::endl << SEPBAR << std::endl << std::endl;

    // std::cout << std::endl << SEPBAR << std::endl << std::endl;

    // // Solve with branch and price
    // std::cout << "Running B&P..." << std::endl;
    // // Copy the original graph
    // Graph gcopy = graph_copy(graph);
    // // Now, execute B&P
    // LP *lp = new LP(gcopy, params, pool, graph, &col2s, true);
    // Node *root = new Node(lp);
    // Col col;
    // BP<Col> bp(col, std::cout, false);
    // bp.set_initial_solution(col2s, col2s.get_n_colors());
    // Stats stats1 = bp.solve(root);
    // stats1.poolSize = pool.size();
    // stats1.print_stats(std::cout);

    // std::cout << std::endl << SEPBAR << std::endl << std::endl;

    // // Solve with compact ilp
    // std::cout << "Running CPLEX with compact ilp formulation..." <<
    // std::endl; Stats stats2 =
    //     solve_ilp(graph, dsaturCol.get_n_colors(), std::cout, dsaturCol);
    // stats2.print_stats(std::cout);

    std::cout << std::endl << SEPBAR << std::endl;
    std::cout << SEPBAR << std::endl;
    std::cout << SEPBAR << std::endl;
    std::cout << SEPBAR << std::endl << std::endl;

    // if (stats1.state == OPTIMAL && stats2.state == OPTIMAL) {
    //   assert(round(stats1.ub) == round(stats2.ub));
    //   assert(round(stats1.lb) == round(stats2.lb));
    // }
  }
}