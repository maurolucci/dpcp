#include "graph.hpp"
#include "lp.hpp"

int main() {

  HyperGraph hg;

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

  ConflictGraph cg;
  init_conflict_graph(hg, cg);

  // Conflict graph:
  // 0: (e1,v1), 1: (e1,v2), 2: (e2,v2), 3: (e2,v3),
  // 4: (e3,v3), 5: (e3,v4), 6: (e4,v4), 7: (e4, v5)
  // 8: (e5,v5), 9: (e5, v1)

  int vcount = 10;
  CGVertex *vlist = new CGVertex[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  int ecount = 0;
  int *elist = new int[2 * num_edges(cg)];
  for (size_t i = 0; i < vcount; ++i) {
    int u = vlist[i];
    for (size_t j = i + 1; j < vcount; ++j) {
      int v = vlist[j];
      auto [e, b] = edge(u, v, cg);
      if (b) {
        elist[2 * ecount] = i;
        elist[2 * ecount + 1] = j;
        ecount++;
      }
    }
  }
  std::cout << "ecount: " << ecount << std::endl;
  std::cout << "elist: ";
  for (int i = 0; i < ecount; ++i) {
    std::cout << elist[2 * i] << " " << elist[2 * i + 1] << " ";
  }
  std::cout << std::endl;

  LP lp(hg, cg, vcount, vlist, ecount, elist, 0);
  lp.optimize();
  std::vector<LP *> branches(2, NULL);
  lp.branch(branches, 2);

  delete branches[0];
  delete branches[1];

  return 0;
}