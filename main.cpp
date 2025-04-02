#include "bp.hpp"
#include "col.hpp"
#include "compact_ilp.hpp"
#include "graph.hpp"
#include "lp.hpp"
#include "stats.hpp"
#include <boost/graph/copy.hpp>

int main() {

  HGraph hg;

  hg.addVertex("v1"); // 0
  hg.addVertex("v2"); // 1
  hg.addVertex("v3"); // 2
  hg.addVertex("v4"); // 3
  hg.addVertex("v5"); // 4

  hg.addHyperedge(hglib::NAME, {"v1", "v2"}); // e1: 0
  hg.addHyperedge(hglib::NAME, {"v2", "v3"}); // e2: 1
  hg.addHyperedge(hglib::NAME, {"v3", "v4"}); // e3: 2
  // Or with ids
  hg.addHyperedge({3, 4}); // e4: 3
  hg.addHyperedge({4, 0}); // e5: 4

  Graph graph;
  get_conflict_graph(hg, graph);

  // Conflict graph:
  // 0: (e1,v1), 1: (e1,v2), 2: (e2,v2), 3: (e2,v3),
  // 4: (e3,v3), 5: (e3,v4), 6: (e4,v4), 7: (e4, v5)
  // 8: (e5,v5), 9: (e5, v1)

  /*
  clear_vertex(1, graph);
  remove_vertex(1, graph);
  clear_vertex(2, graph);
  remove_vertex(2, graph);
  clear_vertex(3, graph);
  remove_vertex(3, graph);
  clear_vertex(4, graph);
  remove_vertex(4, graph);
  clear_vertex(5, graph);
  remove_vertex(5, graph);
  */

  Graph gcopy;
  boost::copy_graph(graph, gcopy);
  Col col(gcopy); // The original graph will be modified during branching
  LP *lp = new LP(graph);
  Node *root = new Node(lp);
  BP<Col> bp(col, std::cout);
  Stats stats1 = bp.solve(root);
  stats1.print_stats(std::cout);

  std::cout << std::endl;
  Stats stats2 = solve_ilp(gcopy, 3, std::cout);
  stats2.print_stats(std::cout);

  return 0;
}