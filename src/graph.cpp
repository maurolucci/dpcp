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
  size_t nVertices = 0;
  for (const auto &e : hg.hyperedges())
    for (const auto &v : *hg.impliedVertices(e->id()))
      add_vertex(VertexInfo{e->id(), v->id(), nVertices++}, graph);
  // Add edge ((e1,v1), (e2,v2)) such that v2 in e1 - v1
  auto [it1, end] = vertices(graph);
  for (; it1 != end; ++it1)
    for (auto it2 = std::next(it1); it2 != end; ++it2) {
      // Unpack PSets
      TypeA e1 = graph[*it1].first, e2 = graph[*it2].first;
      TypeB v1 = graph[*it1].second, v2 = graph[*it2].second;
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

// Vertex branching: v is colored
// Warning: This function modifies the input graph
void vertex_branching1(Graph &graph, Vertex v) {
  TypeA a = graph[v].first;
  TypeA b = graph[v].second;
  auto [it_v, it_end] = vertices(graph);
  for (auto next = it_v; it_v != it_end; it_v = next) {
    next++;
    if (graph[*it_v].first == a && graph[*it_v].second != b) {
      clear_vertex(*it_v, graph);
      remove_vertex(*it_v, graph);
    }
  }
}

// Vertex branching: v is not colored
// Warning: This function modifies the input graph
void vertex_branching2(Graph &graph, Vertex v) {
  clear_vertex(v, graph);
  remove_vertex(v, graph);
}

GraphEnv::GraphEnv(const Graph &graph, Params &params, bool isRoot)
    : GraphEnv(Graph{graph}, params, isRoot){};

GraphEnv::GraphEnv(const Graph &&graph, Params &params, bool isRoot)
    : graph(graph), params(params), getId(), nA(0), nB(0), tyA2idA(), tyB2idB(),
      idA2TyA(), idB2TyB(), snd(), fst(), isRoot(isRoot), isGCP(true),
      isInfeasible(false), isolated() {
  preprocess();
  init_graphenv();
};

GraphEnv::~GraphEnv(){};

void GraphEnv::preprocess() {
  init_preprocess();
  if (params.preprocess[0] == '1' && isRoot)
    preprocess_step1();
  if (params.preprocess[1] == '1')
    preprocess_step2();
  if (params.preprocess[2] == '1')
    preprocess_step3();
}

// Initialize preprocessing
void GraphEnv::init_preprocess() {
  // Initialize tyA2idA, nA (|A|) and snd (V_a, for all a in A)
  // (necessary for preprocessing)
  for (auto v : boost::make_iterator_range(vertices(graph))) {
    TypeA a = graph[v].first;
    bool retA = tyA2idA.insert({a, nA}).second;
    if (retA) {
      snd.push_back(std::vector<Vertex>{v});
      nA++;
    } else
      snd[tyA2idA[a]].push_back(v);
  }
}

// Warring: This function changes the stored graph by adding new edges
void GraphEnv::preprocess_step1() {
  // 1. Make each V_a a clique
  for (size_t id_a = 0; id_a < nA; id_a++) {
    auto it1 = snd[id_a].begin();
    auto it2 = std::next(it1);
    for (; it2 != snd[id_a].end(); ++it1, ++it2)
      if (!edge(*it1, *it2, graph).second)
        add_edge(*it1, *it2, graph);
  }
}

// Warring: This function changes the stored graph by removing edges and
// vertices
void GraphEnv::preprocess_step2() {
  // 2. Iteratively, for all a with |V_a| = 1, wlog V_a = {(a, b_a)},
  // remove N(a,b_a) \cap V_{b_a}. If some V_a = {}, the instance is
  // infeasible.

  // Push candidates (a with |V_a| = 1) into a queue
  std::list<TypeA> queue;
  for (size_t id_a = 0; id_a < nA; id_a++)
    if (snd[id_a].size() == 1)
      queue.push_back(id_a);

  // Process candidates
  while (!queue.empty()) {
    // Pop from queue
    size_t id_a = queue.front();
    queue.pop_front();
    auto v = snd[id_a].front(); // |V_a| = 1
    auto b = graph[v].second;
    // Remove neighbors of v with common b
    auto [it_u, it_end] = adjacent_vertices(v, graph);
    for (auto next = it_u; it_u != it_end; it_u = next) {
      ++next;
      TypeA au = graph[*it_u].first;
      TypeB bu = graph[*it_u].second;
      if (bu == b) {
        // First, remove vertex u from V_{au}
        size_t id_au = tyA2idA[au];
        auto it = std::find(snd[id_au].begin(), snd[id_au].end(), *it_u);
        snd[id_au].erase(it);
        // Check |V_{au}|
        if (snd[id_au].size() == 1)
          // Push u into the queue
          queue.push_back(id_au);
        else if (snd[id_au].size() == 0) {
          // Infeasibility detected
          isInfeasible = true;
          return;
        }
        // Then, remove the vertex from the graph
        clear_vertex(*it_u, graph);
        remove_vertex(*it_u, graph);
      }
    }
  }
}

// Warring: This function changes the stored graph by removing isolated vertices
void GraphEnv::preprocess_step3() {
  // 3. Vanish isolated vertices
  auto [it_u, it_end] = vertices(graph);
  for (auto next = it_u; it_u != it_end; it_u = next) {
    ++next;
    if (degree(*it_u, graph) == 0) {
      auto vi = graph[*it_u];
      isolated.push_back(vi);
      clear_vertex(*it_u, graph);
      remove_vertex(*it_u, graph);
    }
  }
}

void GraphEnv::init_graphenv() {

  // Reset structs intiliazed during preprocessing
  nA = 0;
  tyA2idA.clear();
  snd.clear();

  // Initialize the other struct members
  for (auto v : boost::make_iterator_range(vertices(graph))) {
    getId.insert(std::make_pair(v, getId.size()));
    TypeA a = graph[v].first;
    TypeB b = graph[v].second;
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
}

StableEnv::StableEnv() : stable(), as(), bs(), cost(0.0){};

StableEnv::StableEnv(VertexVector &stable, std::set<TypeA> &as,
                     std::set<TypeB> &bs, double cost)
    : StableEnv(VertexVector{stable}, std::set<TypeA>{as}, std::set<TypeB>{bs},
                cost) {}

StableEnv::StableEnv(VertexVector &&stable, std::set<TypeA> &&as,
                     std::set<TypeB> &&bs, double cost)
    : stable(stable), as(as), bs(bs), cost(cost) {}

bool StableEnv::check(const Graph &graph) {
  for (auto it_v = stable.begin(); it_v != stable.end(); ++it_v)
    for (auto it_u = std::next(it_v); it_u != stable.end(); ++it_u)
      if (edge(*it_v, *it_u, graph).second)
        return false;
  return true;
}