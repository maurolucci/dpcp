#include "graph.hpp"

#include <boost/graph/copy.hpp>

// Graph copying utilities
// Copies the graph structure from src to dst, using vertex2CurrentId to map
// vertex descriptors.
void graph_copy(const Graph& src, const VertexMap<size_t>& vertex2CurrentId,
                Graph& dst) {
  boost::copy_graph(src, dst,
                    boost::vertex_index_map(
                        boost::make_assoc_property_map(vertex2CurrentId)));
}

// Copies the graph structure from src to dst
// Internally, it creates a vertex2CurrentId map that maps vertex descriptors to
// their indices in src.
void graph_copy(const Graph& src, Graph& dst) {
  VertexMap<size_t> vertex2CurrentId;
  for (auto v : boost::make_iterator_range(vertices(src)))
    vertex2CurrentId[v] = vertex2CurrentId.size();
  graph_copy(src, vertex2CurrentId, dst);
}

std::tuple<Graph, Partition, Partition> read_dpcp_instance(
    std::istream& graph, std::istream& partP, std::istream& partQ) {
  Graph g;
  size_t n, m, nP, nQ;
  size_t u, v;
  size_t pi, qj, nPi, nQj;
  char c;
  Partition P, Q;

  graph >> n >> c >> m;

  for (size_t i = 0; i < n; ++i) {
    add_vertex(VertexInfo{i}, g);
  }

  for (size_t i = 0; i < m; ++i) {
    graph >> u >> v;
    add_edge(vertex(u, g), vertex(v, g), g);
  }

  partP >> n >> c >> nP;
  P.resize(nP);
  for (size_t i = 0; i < nP; ++i) {
    partP >> pi >> nPi;
    for (size_t j = 0; j < nPi; ++j) {
      partP >> v;
      P[pi].push_back(vertex(v, g));
    }
  }
  partQ >> n >> c >> nQ;
  Q.resize(nQ);
  for (size_t j = 0; j < nQ; ++j) {
    partQ >> qj >> nQj;
    for (size_t k = 0; k < nQj; ++k) {
      partQ >> v;
      Q[qj].push_back(vertex(v, g));
    }
  }

  return std::make_tuple(std::move(g), std::move(P), std::move(Q));
}

DPCPInst::DPCPInst(const Graph& graph, const Partition& P, const Partition& Q)
    : graph(),
      vertex2CurrentId(),
      P(P),
      Q(Q),
      vertex2Ppart(),
      vertex2Qpart(),
      isolated(),
      isGCP(false),
      isInfeasible(false),
      hasTrivialSolution(false),
      density(0.0) {
  // Copy the graph structure
  graph_copy(graph, this->graph);

  // Initialize vertex2CurrentId, vertex2Ppart and vertex2Qpart maps
  for (Vertex v : boost::make_iterator_range(vertices(this->graph)))
    vertex2CurrentId[v] = vertex2CurrentId.size();
  isGCP = true;
  for (size_t pi = 0; pi < this->P.size(); ++pi) {
    if (P[pi].size() > 1) isGCP = false;
    for (Vertex v : this->P[pi]) vertex2Ppart[v] = pi;
  }
  for (size_t qj = 0; qj < this->Q.size(); ++qj)
    for (Vertex v : this->Q[qj]) vertex2Qpart[v] = qj;

  // Compute density before preprocessing
  size_t n = num_vertices(this->graph);
  if (n > 1) {
    density = static_cast<double>(num_edges(this->graph) * 2) /
              static_cast<double>(n * (n - 1));
  }
}

DPCPInst::DPCPInst(const DPCPInst& dpcp)
    : graph(),
      vertex2CurrentId(),
      P(),
      Q(),
      vertex2Ppart(),
      vertex2Qpart(),
      isGCP(dpcp.is_gcp_instance()),
      isInfeasible(dpcp.is_infeasible_instance()),
      hasTrivialSolution(dpcp.has_trivial_solution()),
      isolated(dpcp.get_isolated_vertices()),
      density(dpcp.get_density()) {
  const Graph& srcGraph = dpcp.get_graph();
  const VertexMap<size_t>& srcVertex2CurrentId = dpcp.get_vertex2CurrentId();
  const Partition& srcP = dpcp.get_P();
  const Partition& srcQ = dpcp.get_Q();

  graph_copy(srcGraph, srcVertex2CurrentId, graph);

  // Map from old vertices (original) to new vertices (copy)
  // Careful: we cannot use vertex descriptors directly, as they may differ
  // between the original and copied graph.
  VertexMap<Vertex> vertexMap;
  for (auto& [vOld, id] : srcVertex2CurrentId) {
    Vertex vNew = vertex(id, graph);
    vertexMap[vOld] = vNew;
    vertex2CurrentId[vNew] = id;
  }

  P.resize(dpcp.get_nP(), std::vector<Vertex>());
  for (size_t pi = 0; pi < dpcp.get_nP(); ++pi)
    for (auto vOld : srcP[pi]) {
      Vertex vNew = vertexMap[vOld];
      P[pi].push_back(vNew);
      vertex2Ppart[vNew] = pi;
    }

  Q.resize(dpcp.get_nQ(), std::vector<Vertex>());
  for (size_t qj = 0; qj < dpcp.get_nQ(); ++qj)
    for (auto vOld : srcQ[qj]) {
      Vertex vNew = vertexMap[vOld];
      Q[qj].push_back(vNew);
      vertex2Qpart[vNew] = qj;
    }
}

DPCPInst::~DPCPInst() {}

GCPGraph DPCPInst::get_gcp_graph() const {
  GCPGraph gcpGraph;
  // Add vertices
  for (size_t qj = 0; qj < get_nQ(); ++qj) add_vertex(qj, gcpGraph);
  // Add edges
  for (auto e : boost::make_iterator_range(edges(graph))) {
    auto b1 = get_Q_part(source(e, graph));
    auto b2 = get_Q_part(target(e, graph));
    if (!edge(b1, b2, gcpGraph).second) add_edge(b1, b2, gcpGraph);
  }
  return gcpGraph;
}

void DPCPInst::remove_vertex(Vertex v) {
  size_t pi = get_P_part(v);
  size_t qj = get_Q_part(v);

  // First, we remove v from its P-part and Q-part. If any of the parts becomes
  // empty, we need to remove it and update the corresponding maps as well.

  auto it = std::find(P[pi].begin(), P[pi].end(), v);
  assert(it != P[pi].end());  // v should be in P[pi]
  P[pi].erase(it);
  if (P[pi].empty()) {
    isInfeasible = true;
    for (auto it_i = pi + 1; it_i < P.size(); ++it_i) {
      P[it_i - 1] = std::move(P[it_i]);
      for (auto v : P[it_i - 1]) vertex2Ppart[v] = it_i - 1;
    }
    P.pop_back();
  }

  it = std::find(Q[qj].begin(), Q[qj].end(), v);
  assert(it != Q[qj].end());  // v should be in Q[qj]
  Q[qj].erase(it);
  if (Q[qj].empty()) {
    for (auto it_j = qj + 1; it_j < Q.size(); ++it_j) {
      Q[it_j - 1] = std::move(Q[it_j]);
      for (auto v : Q[it_j - 1]) vertex2Qpart[v] = it_j - 1;
    }
    Q.pop_back();
  }

  vertex2Ppart.erase(v);
  vertex2Qpart.erase(v);

  // Finally, we remove the vertex from the graph and update vertex2CurrentId

  clear_vertex(v, graph);
  boost::remove_vertex(v, graph);

  vertex2CurrentId.clear();
  for (auto v : boost::make_iterator_range(vertices(graph)))
    vertex2CurrentId[v] = vertex2CurrentId.size();
}

void DPCPInst::preprocess(bool clique) {
  if (clique) preprocess_step1();
  preprocess_step2();
  preprocess_step3();
  preprocess_step4();
}

// Preprocess #1: Make each P-part a clique by adding missing edges.
void DPCPInst::preprocess_step1() {
  isGCP = true;
  for (size_t pi = 0; pi < get_nP(); ++pi) {
    for (auto it_u = P[pi].begin(); it_u != P[pi].end(); ++it_u)
      for (auto it_v = std::next(it_u); it_v != P[pi].end(); ++it_v)
        if (!edge(*it_u, *it_v, graph).second) add_edge(*it_u, *it_v, graph);
  }
}

// Preprocess #2: If a P-part has only one vertex v, we remove their neighbors
// in the same Q-part as v.
void DPCPInst::preprocess_step2() {
  for (size_t pi = 0; pi < get_nP(); ++pi)
    if (P[pi].size() == 1) {
      Vertex v = P[pi].front();
      bool removed = false;
      auto [it_u, it_end] = adjacent_vertices(v, graph);
      for (auto next = it_u; it_u != it_end; it_u = next) {
        ++next;
        Vertex u = *it_u;
        if (get_Q_part(u) == get_Q_part(v)) {
          remove_vertex(u);
          removed = true;
        }
      }
      // The removal of u may create new P-parts of size 1, so we need to repeat
      // the process until no more vertices can be removed.
      if (removed) return preprocess_step2();
    }
}

// Preprocess #3: Remove isolated vertices.
void DPCPInst::preprocess_step3() {
  bool oldFeasibility = isInfeasible;
  auto [it_u, it_end] = vertices(graph);
  for (auto next = it_u; it_u != it_end; it_u = next) {
    ++next;
    Vertex u = *it_u;
    if (degree(u, graph) == 0) {
      // Since we complete cliques, u should be the only vertex in its P-part
      assert(P[get_P_part(u)].size() == 1);
      // Save the original id of u
      isolated.push_back(IsolatedVertex{graph[u].id});
      remove_vertex(u);
    }
  }
  // The removal of a vertex can raise infeasibility, but the removal of
  // isolated vertices do not change the feasibility of the instance, so we
  // restore the original feasibility value.
  isInfeasible = oldFeasibility;
}

// Preprocess #4: If there is only one P-part, then we have a trivial
// solution.
void DPCPInst::preprocess_step4() {
  if (get_nP() == 1) hasTrivialSolution = true;
}

void DPCPInst::preselect_vertex(Vertex v) {
  size_t pi = get_P_part(v);
  assert(P[pi].size() > 1);
  // Local copy of the P-part, as we will modify it during the loop
  VertexVector Pi(P[pi]);
  for (auto u : Pi)
    if (u != v) remove_vertex(u);
}

void DPCPInst::forbid_vertex(Vertex v) { remove_vertex(v); }

StableEnv::StableEnv() : stable(), ps(), qs(), cost(0.0) {}

StableEnv::StableEnv(VertexVector& stable, std::set<size_t>& ps,
                     std::set<size_t>& qs, double cost)
    : StableEnv(VertexVector{stable}, std::set<size_t>{ps},
                std::set<size_t>{qs}, cost) {}

StableEnv::StableEnv(VertexVector&& stable, std::set<size_t>&& ps,
                     std::set<size_t>&& qs, double cost)
    : stable(std::move(stable)),
      ps(std::move(ps)),
      qs(std::move(qs)),
      cost(cost) {}

bool StableEnv::check(const Graph& graph) {
  for (auto it_v = stable.begin(); it_v != stable.end(); ++it_v)
    for (auto it_u = std::next(it_v); it_u != stable.end(); ++it_u)
      if (edge(*it_v, *it_u, graph).second) return false;
  return true;
}

void StableEnv::add_vertex(const Vertex v, size_t pi, size_t qj) {
  stable.push_back(v);
  ps.insert(pi);
  qs.insert(qj);
}
