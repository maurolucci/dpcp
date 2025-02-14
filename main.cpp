#include "graph.hpp"
#include "lp.hpp"

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

  std::cout << num_vertices(graph) << " " << num_edges(graph) << std::endl;

  LP lp(graph);
  LP_STATE state = lp.optimize();
  (void)state;

  std::vector<LP *> branches(2, NULL);
  lp.branch(branches, 5);
  // delete branches[0];
  // delete branches[1];

  return 0;
}