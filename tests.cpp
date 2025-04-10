#include "bp.hpp"
#include "col.hpp"
#include "compact_ilp.hpp"
#include "graph.hpp"
#include "heur.hpp"
#include "lp.hpp"
#include "random.hpp"
#include "stats.hpp"

#include <algorithm> // shuffle
#include <boost/graph/copy.hpp>
#include <filesystem>
#include <iostream>
#include <random> // std::default_random_engine

using recursive_directory_iterator =
    std::filesystem::recursive_directory_iterator;

std::string SEPBAR =
    "*************************************************************";

int main() {

  // Set seed
  set_seed(0);

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

    // Build conflict graph
    Graph graph;
    get_conflict_graph(hg, graph);
    std::cout << "Conflict graph:" << std::endl;
    std::cout << "\tVertices: " << num_vertices(graph) << std::endl;
    std::cout << "\tEdges: " << num_edges(graph) << std::endl;

    std::cout << std::endl << SEPBAR << std::endl << std::endl;

    // Copy conflict graph
    Graph gcopy;
    boost::copy_graph(graph, gcopy);

    // Solve with heuristic
    GraphEnv genv(gcopy);
    Col dsaturCol;
    std::cout << "Running heuristic..." << std::endl;
    Stats stats0 = heur_solve(genv, genv.idA2TyA, dsaturCol, 200);
    stats0.print_stats(std::cout);

    std::cout << std::endl << SEPBAR << std::endl << std::endl;

    // Solve with compact ilp
    std::cout << "Running CPLEX with compact ilp formulation..." << std::endl;
    Stats stats2 =
        solve_ilp(gcopy, dsaturCol.get_n_colors(), std::cout, dsaturCol);
    stats2.print_stats(std::cout);

    std::cout << std::endl << SEPBAR << std::endl << std::endl;

    // Solve with branch and price
    std::cout << "Running B&P..." << std::endl;
    Col col;
    LP *lp = new LP(graph);
    Node *root = new Node(lp);
    BP<Col> bp(col, std::cout, false);
    Stats stats1 = bp.solve(root);
    stats1.print_stats(std::cout);

    std::cout << std::endl << SEPBAR << std::endl << std::endl;

    if (stats1.state == OPTIMAL && stats2.state == OPTIMAL) {
      assert(round(stats1.ub) == round(stats2.ub));
      assert(round(stats1.lb) == round(stats2.lb));
    }
  }
}