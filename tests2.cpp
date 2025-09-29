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

  for (const auto &file : recursive_directory_iterator("input/infeasibles/")) {

    if (file.path().extension() != ".graph")
      continue;

    std::cout << "*** Solving: " << file.path().string() << " ***" << std::endl;

    // Open files
    std::ifstream inGraph(file.path().string());
    std::ifstream inPartA(file.path().parent_path().string() + "/" +
                          file.path().stem().string() + ".partA");
    std::ifstream inPartB(file.path().parent_path().string() + "/" +
                          file.path().stem().string() + ".partB");
    if (!inGraph.is_open() || !inPartA.is_open() || !inPartB.is_open()) {
      std::cerr << "Error opening files" << std::endl;
      return -1;
    }

    // Read DPCP instance
    Graph dpcp = read_dpcp_instance(inGraph, inPartA, inPartB);
    std::cout << "DPCP instance:" << std::endl;
    std::cout << "\tVertices: " << num_vertices(dpcp) << std::endl;
    std::cout << "\tEdges: " << num_edges(dpcp) << std::endl;

    // Read parameters
    Params params;

    std::cout << std::endl << SEPBAR << std::endl << std::endl;

    // Solve with branch and price
    std::cout << "Running B&P..." << std::endl;
    // Copy the original graph
    Graph gcopy = graph_copy(dpcp);
    // Now, execute B&P
    Pool pool;
    LP *lp = new LP(gcopy, params, pool, dpcp, true);
    Node *root = new Node(lp);
    Col col;
    BP<Col> bp(params, std::cout, col);
    Stats stats1 = bp.solve(root);
    stats1.poolSize = pool.size();
    stats1.print_stats(std::cout);

    std::cout << std::endl << SEPBAR << std::endl << std::endl;

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