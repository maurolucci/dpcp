#include "graph.hpp"

#include <boost/graph/copy.hpp>

// Read a DPCP instance from input stream
// Returns the graph and the sizes of the two partitions
std::tuple<Graph, size_t, size_t> read_dpcp_instance(std::istream &graph,
                                                     std::istream &partA,
                                                     std::istream &partB) {
  Graph g;
  size_t n, m, nA, nB;
  size_t a, nVa, nVb;
  char c;

  // Read the partitions
  std::map<size_t, size_t> vertexToA, vertexToB;
  partA >> n >> c >> nA;
  for (size_t i = 0; i < nA; ++i) {
    partA >> a >> nVa;
    for (size_t j = 0; j < nVa; ++j) {
      size_t v;
      partA >> v;
      vertexToA[v] = a;
    }
  }
  partB >> n >> c >> nB;
  for (size_t i = 0; i < nB; ++i) {
    partB >> a >> nVb;
    for (size_t j = 0; j < nVb; ++j) {
      size_t v;
      partB >> v;
      vertexToB[v] = a;
    }
  }

  // Create vertices with their partition info
  std::vector<Vertex> vertices(n);
  for (size_t i = 0; i < n; ++i) {
    TypeA a = vertexToA[i];
    TypeB b = vertexToB[i];
    vertices[i] = add_vertex(VertexInfo{a, b, i}, g);
  }

  // Add edges
  graph >> n >> c >> m;
  for (size_t i = 0; i < m; ++i) {
    size_t u, v;
    graph >> u >> v;
    add_edge(vertices[u], vertices[v], g);
  }

  return std::make_tuple(g, nA, nB);
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

void graph_copy(const Graph &src, const std::map<Vertex, size_t> &getId,
                Graph &dst) {
  boost::copy_graph(
      src, dst, boost::vertex_index_map(boost::make_assoc_property_map(getId)));
}

void graph_copy(const Graph &src, Graph &dst) {
  // boost::copy_graph needs a map from Vertex to size_t
  std::map<Vertex, size_t> index;
  for (auto v : boost::make_iterator_range(vertices(src)))
    index.insert(std::make_pair(v, index.size()));
  graph_copy(src, index, dst);
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

GraphEnv::GraphEnv(Graph *graph, bool preprocess1, bool preprocess2,
                   bool preprocess3, bool preprocess4, bool isRoot)
    : graphPtr(graph), graph(*graph), getId(), nA(0), nB(0), Va(), Vb(),
      tyA2idA(), tyB2idB(), idA2TyA(), idB2TyB(), snd(), isRoot(isRoot),
      isGCP(true), isInfeasible(false), hasTrivialSolution(false), isolated() {
  // Density
  density = static_cast<double>((num_edges(*graph)) * 2) /
            (num_vertices(*graph) * (num_vertices(*graph) - 1));
  // First, apply preprocessing
  preprocess(preprocess1, preprocess2, preprocess3, preprocess4);
  // Then, initialize the other struct members
  init_graphenv();
};

GraphEnv::~GraphEnv(){};

void GraphEnv::preprocess(bool preprocess1, bool preprocess2, bool preprocess3,
                          bool preprocess4) {
  init_preprocess();
  if (preprocess1 && isRoot)
    preprocess_step1();
  if (preprocess2)
    preprocess_step2();
  if (preprocess3)
    preprocess_step3();
  if (preprocess4)
    preprocess_step4();
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
        auto vu = *it_u;
        clear_vertex(vu, graph);
        remove_vertex(vu, graph);
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
      nA--;
    }
  }
}

void GraphEnv::preprocess_step4() {
  // 4. Check if n = 1
  if (nA == 1)
    hasTrivialSolution = true;
}

void GraphEnv::init_graphenv() {

  // Reset structs intiliazed during preprocessing
  nA = 0;
  tyA2idA.clear();
  snd.clear();

  std::set<TypeB> B;

  // Initialize the other struct members
  for (auto v : boost::make_iterator_range(vertices(graph))) {
    getId.insert(std::make_pair(v, getId.size()));
    TypeA a = graph[v].first;
    TypeB b = graph[v].second;
    bool retA = tyA2idA.insert({a, nA}).second;
    // bool retB = tyB2idB.insert({b, nB}).second;
    bool retB = B.insert(b).second;
    if (retA) {
      idA2TyA.push_back(a);
      snd.push_back(std::vector<Vertex>{v});
      Va.emplace(a, std::vector<Vertex>{v});
      nA++;
    } else {
      snd[tyA2idA[a]].push_back(v);
      Va[a].push_back(v);
      isGCP = false;
    }
    if (retB) {
      Vb.emplace(b, std::vector<Vertex>{v});
      nB++;
    } else {
      Vb[b].push_back(v);
    }
  }

  // Use B to initializa idB2TyB
  for (const auto &b : B)
    idB2TyB.push_back(b);

  // Sort idB2TyB by decreasing order of |Vb|
  std::sort(idB2TyB.begin(), idB2TyB.end(),
            [this](const TypeB b1, const TypeB b2) {
              return Vb[b1].size() > Vb[b2].size();
            });

  // Initialize tyB2idB
  for (size_t i = 0; i < idB2TyB.size(); ++i)
    tyB2idB[idB2TyB[i]] = i;
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

void StableEnv::add_vertex(const Vertex v, const TypeA a, const TypeB b) {
  stable.push_back(v);
  as.insert(a);
  bs.insert(b);
}
