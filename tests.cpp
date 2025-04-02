#include "bp.hpp"
#include "col.hpp"
#include "compact_ilp.hpp"
#include "graph.hpp"
#include "lp.hpp"
#include "stats.hpp"

#include <boost/graph/copy.hpp>
#include <filesystem>
#include <iostream>

using recursive_directory_iterator =
    std::filesystem::recursive_directory_iterator;

int main() {

  for (const auto &file : recursive_directory_iterator("input/tests/")) {

    if (file.path().extension() != ".txt")
      continue;

    std::cout << "\n\n\n*** Solving: " << file << "***" << std::endl;

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

    // Copy conflict graph
    Graph gcopy;
    boost::copy_graph(graph, gcopy);

    // Solve with branch and price
    Col col(gcopy);
    LP *lp = new LP(graph);
    Node *root = new Node(lp);
    BP<Col> bp(col, std::cout);
    Stats stats1 = bp.solve(root);
    stats1.print_stats(std::cout);

    // Solve with compact ilp
    std::cout << std::endl;
    Stats stats2 = solve_ilp(gcopy, 5, std::cout);
    stats2.print_stats(std::cout);

    if (stats1.state == OPTIMAL && stats2.state == OPTIMAL) {
      assert(stats1.ub == stats2.ub);
      assert(stats1.lb == stats2.lb);
    }
  }
}