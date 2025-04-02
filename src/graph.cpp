#include "graph.hpp"

void read_hypergrah(HGraph &hg, std::istream &input) {

  size_t n, m, mm, v;
  input >> n >> m;

  // Add vertices
  for (size_t v = 0; v < n; ++v)
    hg.addVertex();

  // Add hyperedges
  for (size_t e = 0; e < m; ++e) {
    input >> mm;
    std::vector<HVertex> vertices;
    for (size_t i = 0; i < mm; ++i) {
      input >> v;
      vertices.push_back(v);
    }
    hg.addHyperedge(vertices);
  }
}

void get_conflict_graph(const HGraph &hg, Graph &graph) {
  // Add vertices (pointed sets)
  for (const auto &e : hg.hyperedges())
    for (const auto &v : *hg.impliedVertices(e->id()))
      add_vertex(std::make_pair(e->id(), v->id()), graph);
  // Add edge ((e1,v1), (e2,v2)) such that v2 in e1 - v1
  auto [it1, end] = vertices(graph);
  for (; it1 != end; ++it1)
    for (auto it2 = std::next(it1); it2 != end; ++it2) {
      // Unpack PSets
      auto [e1, v1] = graph[*it1];
      auto [e2, v2] = graph[*it2];
      if (v1 != v2 &&
          (hg.isVertexOfHyperedge(v2, e1) || hg.isVertexOfHyperedge(v1, e2)))
        add_edge(*it1, *it2, graph);
    }
  return;
}

void get_gcp_graph(Graph &src, GCPGraph &dst, std::map<TypeB, size_t> &tyB2idB,
                   std::vector<TypeB> &idB2TyB) {
  // Add vertices
  for (auto b : idB2TyB)
    add_vertex(b, dst);
  // Add edges
  for (auto e : boost::make_iterator_range(edges(src))) {
    auto b1 = src[source(e, src)].second;
    auto b2 = src[target(e, src)].second;
    if (!edge(tyB2idB[b1], tyB2idB[b2], dst).second)
      add_edge(tyB2idB[b1], tyB2idB[b2], dst);
  }
}