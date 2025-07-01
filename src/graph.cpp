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

GraphEnv::GraphEnv(const Graph &graph) : GraphEnv(Graph{graph}){};

GraphEnv::GraphEnv(const Graph &&graph) : graph(graph) {

  // Fill maps
  nA = 0, nB = 0;
  isGCP = true;
  for (auto v : boost::make_iterator_range(vertices(graph))) {
    auto [a, b] = graph[v];
    bool retA = tyA2idA.insert({a, nA}).second;
    bool retB = tyB2idB.insert({b, nB}).second;
    if (retA) {
      idA2TyA.push_back(a);
      snd.push_back(std::vector<Vertex>{v});
      nA++;
    } else {
      snd[tyA2idA[a]].push_back(v);
      isGCP = false;
    }
    if (retB) {
      idB2TyB.push_back(b);
      fst.push_back(std::vector<Vertex>{v});
      nB++;
    } else
      fst[tyB2idB[b]].push_back(v);
  }
};

GraphEnv::~GraphEnv(){};

StableEnv::StableEnv() : stable(), as(), bs(), cost(0.0){};

StableEnv::StableEnv(VertexVector &stable, std::set<TypeB> &as,
                     std::set<TypeB> &bs, double cost)
    : StableEnv(VertexVector{stable}, VertexSet{as}, VertexSet{bs}, cost) {}

StableEnv::StableEnv(VertexVector &&stable, std::set<TypeB> &&as,
                     std::set<TypeB> &&bs, double cost)
    : stable(stable), as(as), bs(bs), cost(cost) {}

bool StableEnv::check(const Graph &graph) {
  for (auto it_v = stable.begin(); it_v != stable.end(); ++it_v)
    for (auto it_u = std::next(it_v); it_u != stable.end(); ++it_u)
      if (edge(*it_v, *it_u, graph).second)
        return false;
  return true;
}