#include "graph.hpp"

void init_conflict_graph(const HyperGraph &hg, ConflictGraph &cg) {
  // Add vertices (PSet)
  for (const auto &e : hg.hyperedges())
    for (const auto &v : *hg.impliedVertices(e->id()))
      add_vertex(std::make_pair(e->id(), v->id()), cg);
  // Add edge ((e1,v1), (e2,v2)) such that v2 in e1 - v1
  auto [it1, end] = vertices(cg);
  for (; it1 != end; ++it1)
    for (auto it2 = std::next(it1); it2 != end; ++it2) {
      // Unpack PSets
      auto [e1, v1] = cg[*it1];
      auto [e2, v2] = cg[*it2];
      if (v1 != v2 &&
          (hg.isVertexOfHyperedge(v2, e1) || hg.isVertexOfHyperedge(v1, e2)))
        add_edge(*it1, *it2, cg);
    }
  return;
}