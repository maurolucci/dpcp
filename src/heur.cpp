#include "heur.hpp"

#include "graph.hpp"
#include "random.hpp"
extern "C" {
#include "color.h"
}

#include <chrono>
#include <fstream>
#include <limits>
#include <queue>
#include <random>
#include <ranges>
#include <stdexcept>

using ClockType = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::_V2::system_clock::time_point;

/******************* */
/* 2-STEP HEURISTICS */
/******************* */

// Two-step greedy heuristic criteria

// Criteria 0: DEG-REAL
// Degree of v in G, ignoring neighbors in P_{i(v)} and Q_{j(v)}
size_t get_real_degree(const DPCPInst& dpcp,
                       const std::map<Vertex, bool>& removed,
                       const VertexVector& selected,
                       const std::map<size_t, std::set<size_t>>& adj,
                       const Vertex& v) {
  size_t degree = 0;
  for (Vertex u :
       boost::make_iterator_range(adjacent_vertices(v, dpcp.get_graph()))) {
    if (removed.at(u)) continue;
    if (dpcp.get_P_part(u) != dpcp.get_P_part(v) &&
        dpcp.get_Q_part(u) != dpcp.get_Q_part(v))
      degree++;
  }
  return degree;
}

// Criteria 1: DEG-Q (used for breaking ties)
// Degree of v in G[Q_{j(v)}]
size_t get_Q_degree(const DPCPInst& dpcp, const std::map<Vertex, bool>& removed,
                    const VertexVector& selected,
                    const std::map<size_t, std::set<size_t>>& adj,
                    const Vertex& v) {
  size_t neighbors = 0;
  for (Vertex u :
       boost::make_iterator_range(adjacent_vertices(v, dpcp.get_graph()))) {
    if (removed.at(u)) continue;
    if (dpcp.get_Q_part(u) == dpcp.get_Q_part(v)) neighbors++;
  }
  return neighbors;
}

// Criteria 2: DEG-COLLAPSED
// Degree of v in the collapsed graph according to Q,
// ignoring neighbors in P_{i(v)} and Q_{j(v)}
size_t get_collapsed_degree(const DPCPInst& dpcp,
                            const std::map<Vertex, bool>& removed,
                            const VertexVector& selected,
                            const std::map<size_t, std::set<size_t>>& adj,
                            const Vertex& v) {
  std::set<size_t> Qdegree;
  for (Vertex u :
       boost::make_iterator_range(adjacent_vertices(v, dpcp.get_graph()))) {
    if (removed.at(u)) continue;
    if (dpcp.get_P_part(u) != dpcp.get_P_part(v) &&
        dpcp.get_Q_part(u) != dpcp.get_Q_part(v))
      Qdegree.insert(dpcp.get_Q_part(u));
  }
  return Qdegree.size();
}

// Criteria 3: EDG
// No. of edges added to \tilde{G}[W]
size_t get_n_new_edge(const DPCPInst& dpcp,
                      const std::map<Vertex, bool>& removed,
                      const VertexVector& selected,
                      const std::map<size_t, std::set<size_t>>& adj,
                      const Vertex& v) {
  std::set<size_t> adjBs;
  size_t bv = dpcp.get_Q_part(v);
  for (Vertex u : selected) {
    if (!edge(v, u, dpcp.get_graph()).second) continue;
    size_t bu = dpcp.get_Q_part(u);
    if (bv == bu) continue;
    if (bv < bu && (!adj.contains(bv) || !adj.at(bv).contains(bu)))
      adjBs.insert(bu);
    else if (bv > bu && (!adj.contains(bu) || !adj.at(bu).contains(bv)))
      adjBs.insert(bu);
  }
  return adjBs.size();
}

size_t evaluate_vertex(const DPCPInst& dpcp,
                       const std::map<Vertex, bool>& removed,
                       const VertexVector& selected,
                       const std::map<size_t, std::set<size_t>>& adj,
                       size_t variant, const Vertex& v) {
  if (variant == 0)
    return get_real_degree(dpcp, removed, selected, adj, v);
  else if (variant == 1)
    return get_Q_degree(dpcp, removed, selected, adj, v);
  else if (variant == 2)
    return get_collapsed_degree(dpcp, removed, selected, adj, v);
  else if (variant == 3)
    return get_n_new_edge(dpcp, removed, selected, adj, v);
  else {
    if (dpcp.get_density() <= 0.6)
      return get_n_new_edge(dpcp, removed, selected, adj, v);
    else
      return get_collapsed_degree(dpcp, removed, selected, adj, v);
  }
}

// Function to select the next vertex v following a greedy strategy
Vertex greedy_vertex_selector(const DPCPInst& dpcp,
                              const VertexVector& candidates,
                              const std::map<Vertex, bool>& removed,
                              const VertexVector& selected,
                              std::map<size_t, std::set<size_t>>& adj,
                              const Params& params) {
  Vertex minVertex = NULL;
  size_t minVal = std::numeric_limits<size_t>::max();
  size_t minTie = std::numeric_limits<size_t>::max();
  for (Vertex v : candidates) {
    if (removed.at(v)) continue;
    size_t val = evaluate_vertex(dpcp, removed, selected, adj,
                                 params.heuristic2stepVariant, v);
    if (val < minVal) {
      minVal = val;
      minTie = evaluate_vertex(dpcp, removed, selected, adj, 1, v);
      minVertex = v;
    } else if (val == minVal) {
      size_t tie = evaluate_vertex(dpcp, removed, selected, adj, 1, v);
      if (tie < minTie) {
        minTie = tie;
        minVertex = v;
      }
    }
  }
  return minVertex;
}

// Function to select the next vertex v following a semigreedy strategy
// We use a restricted candidate list (RCL) of size |RCL| = max(1, alpha*|C|)
// where C is the set of candidates and alpha in [0,1] is a parameter
// We store in the RCL the best candidates and we select randomly one of them
Vertex semigreedy_vertex_selector(const DPCPInst& dpcp,
                                  const VertexVector& candidates,
                                  const std::map<Vertex, bool>& removed,
                                  const VertexVector& selected,
                                  std::map<size_t, std::set<size_t>>& adj,
                                  const Params& params) {
  dpcp.check_consistency();
  // First, find the lowest and highest degree among the candidates
  size_t minVal = std::numeric_limits<size_t>::max();
  size_t maxVal = 0;
  std::map<Vertex, size_t> valMap;

  // Print the candidates and their values
  std::cout << "Candidates: " << std::endl;
  for (Vertex v : candidates) std::cout << v << std::endl;

  for (Vertex v : candidates) {
    if (removed.at(v)) continue;
    size_t val = evaluate_vertex(dpcp, removed, selected, adj,
                                 params.heuristic2stepVariant, v);
    valMap[v] = val;
    if (val < minVal) minVal = val;
    if (val > maxVal) maxVal = val;
  }
  double cutoff = minVal + params.heuristicSemigreedyAlpha * (maxVal - minVal);

  // Build the RCL
  VertexVector rcl;
  for (auto [v, d] : valMap)
    if (d <= cutoff) rcl.push_back(v);

  if (rcl.empty()) return NULL;

  // Select randomly a vertex from the RCL
  auto r = rand_int(rng) % rcl.size();

  return rcl[r];
}

using Heur2SVertexSelector = Vertex (*)(const DPCPInst&, const VertexVector&,
                                        const std::map<Vertex, bool>&,
                                        const VertexVector&,
                                        std::map<size_t, std::set<size_t>>&,
                                        const Params& params);

// Update the information after removing a vertex v
void update_info(const DPCPInst& dpcp, std::map<Vertex, bool>& removed,
                 Vertex v, std::map<size_t, size_t>& nP) {
  removed.at(v) = true;
  size_t pi = dpcp.get_P_part(v);
  nP.at(pi)--;
}

// First step of the 2-step heuristic for DPCP
// Select a set of vertices of size |P|, one from each P[pi], such that the
// selected vertices of each Q[qj] are a stable set
bool first_step(const DPCPInst& dpcp, VertexVector& selected,
                std::map<size_t, std::set<size_t>>& adj, const Params& params,
                Heur2SVertexSelector vertexSelector) {
  // Print the candidates and their values
  std::cout << "Vertices: " << std::endl;
  for (Vertex v : boost::make_iterator_range(vertices(dpcp.get_graph())))
    std::cout << v << std::endl;

  // Map with the removed vertices
  std::map<Vertex, bool> removed;
  for (Vertex u : boost::make_iterator_range(vertices(dpcp.get_graph())))
    removed.emplace(u, false);

  // Fill the map with the size of each P[pi]
  std::map<size_t, size_t> nP;
  for (size_t pi = 0; pi < dpcp.get_nP(); ++pi)
    nP.emplace(pi, dpcp.get_P()[pi].size());

  // Subset of P indices that has not been processed yet
  std::set<size_t> unprocessedP;
  for (size_t pi = 0; pi < dpcp.get_nP(); ++pi) unprocessedP.insert(pi);

  // Comparison function for selecting the next pi to process
  std::function<bool(const size_t, const size_t)> compareFunc;
  if (vertexSelector == greedy_vertex_selector) {
    // For greedy strategy, choose pi such that P[pi] has currently the
    // minimum size, break ties by pi index
    compareFunc = [&nP](const size_t pi1, const size_t pi2) {
      return (nP.at(pi1) < nP.at(pi2) ||
              (nP.at(pi1) == nP.at(pi2) && pi1 < pi2));
    };
  } else {
    // Otherwise, choose pi such that P[pi] has currently the minimum size,
    // but break ties randomly
    compareFunc = [&nP](const size_t pi1, const size_t pi2) {
      return (nP.at(pi1) < nP.at(pi2) ||
              (nP.at(pi1) == nP.at(pi2) && (rand_int(rng) % 2) == 0));
    };
  }

  while (selected.size() < dpcp.get_nP()) {
    // Choose an unprocessed pi
    size_t pi = *std::min_element(unprocessedP.begin(), unprocessedP.end(),
                                  compareFunc);
    unprocessedP.erase(pi);

    // Choose a vertex, with some criterion
    Vertex v =
        vertexSelector(dpcp, dpcp.get_P()[pi], removed, selected, adj, params);
    // info.at(v).print_info();
    if (v == NULL) return false;  // No vertex can be selected from P[pi]

    // Get label qj of v
    size_t qj = dpcp.get_Q_part(v);

    // Update adjacent list and selected
    for (Vertex u : selected) {
      if (!edge(u, v, dpcp.get_graph()).second) continue;
      size_t qj2 = dpcp.get_Q_part(u);
      if (qj < qj2)
        adj[qj].insert(qj2);
      else if (qj > qj2)
        adj[qj2].insert(qj);
    }
    selected.push_back(v);

    // Remove the vertices invalidated by choosing v, i.e.
    // vertices of P[pi] \ v and neighbors of v that also belong to Q[qj]
    for (Vertex u : dpcp.get_P()[pi])
      if (!removed.at(u) && u != v) update_info(dpcp, removed, u, nP);
    for (Vertex u :
         boost::make_iterator_range(adjacent_vertices(v, dpcp.get_graph())))
      if (!removed.at(u) && dpcp.get_Q_part(u) == qj)
        update_info(dpcp, removed, u, nP);
  }

  return selected.size() == dpcp.get_nP();
}

// Second step of the 2-step heuristic for DPCP
// Color the selected vertices using DSATUR on the graph induced by them,
// collapsing vertices of each Q[qj]
void dpcp_dsatur_heur(const DPCPInst& dpcp, VertexVector& selected,
                      std::map<size_t, std::set<size_t>>& adj, Col& col) {
  std::map<size_t, size_t>
      qjToSubIndex;  // Map from qj to its new index (in the subgraph)
  std::vector<size_t>
      subIndexToQj;  // Map from new index (in the subgraph) to qj
  std::map<size_t, VertexVector>
      repr;        // Map from qj to the vector of vertices that qj represents
  int ecount = 0;  // Number of edges in the subgraph
  int elist[2 * num_edges(dpcp.get_graph())];  // Edgle list of the subgraph

  // First, initialize the qjToSubIndex, subIndexToQj and repr structures
  // This will be used for building the vertex set of the subgraph to color
  for (Vertex v : selected) {
    size_t qj1 = dpcp.get_Q_part(v);
    if (!qjToSubIndex.contains(qj1)) {
      qjToSubIndex[qj1] = qjToSubIndex.size();
      subIndexToQj.push_back(qj1);
      repr[qj1] = std::vector<Vertex>();
    }
    repr[qj1].push_back(v);
  }

  // Seconds, build the edge list of the subgraph to color
  for (auto& [qj1, adjQj] : adj)
    for (size_t qj2 : adjQj) {
      elist[2 * ecount] = qjToSubIndex[qj1];
      elist[2 * ecount + 1] = qjToSubIndex[qj2];
      ecount++;
    }

  // Third, color the subgraph using DSATUR
  int ncolors = 0;
  COLORset* colorclasses = NULL;
  COLORdsatur(qjToSubIndex.size(), ecount, elist, &ncolors, &colorclasses);

  // Recover coloring
  col.reset_coloring();
  for (int k = 0; k < ncolors; ++k)
    for (int j = 0; j < colorclasses[k].count; ++j) {
      size_t qj = subIndexToQj[colorclasses[k].members[j]];
      for (Vertex v : repr[qj]) col.set_color(dpcp, dpcp.get_current_id(v), k);
    }
}

// General two-step greedy heuristic for DPCP
HeurStats dpcp_2_step_greedy_heur(const DPCPInst& dpcp, Col& col,
                                  const Params& params) {
  TimePoint start = ClockType::now();
  HeurStats stats;

  // First step
  VertexVector selected;                   // Vector of selected vertices
  std::map<size_t, std::set<size_t>> adj;  // Adjacent list of the subgraph
                                           // induced by the selected vertices
  bool success =
      first_step(dpcp, selected, adj, params, greedy_vertex_selector);
  if (!success) {
    stats.state = UNKNOWN;
  } else {
    // Second step
    dpcp_dsatur_heur(dpcp, selected, adj, col);
    assert(col.check_coloring(dpcp));
    stats.state = FEASIBLE;
    stats.value = static_cast<double>(col.get_n_colors());
    stats.bestTime =
        std::chrono::duration<double>(ClockType::now() - start).count();
    stats.bestIter = 0;
  }

  TimePoint end = ClockType::now();
  stats.totalTime = std::chrono::duration<double>(end - start).count();
  stats.totalIters = 1;

  return stats;
}

// General two-step semigreedy heuristic for DPCP
HeurStats dpcp_2_step_semigreedy_heur(const DPCPInst& dpcp, Col& col,
                                      const Params& params,
                                      std::ostream& iterFile) {
  dpcp.check_consistency();
  TimePoint start = ClockType::now();
  HeurStats stats;
  stats.totalIters =
      params.heuristicSemigreedyIter * num_vertices(dpcp.get_graph());

  for (size_t i = 0; i < stats.totalIters; ++i) {
    // First step
    VertexVector selected;                   // Vector of selected vertices
    std::map<size_t, std::set<size_t>> adj;  // Adjacent list of the subgraph
                                             // induced by the selected vertices
    bool success =
        first_step(dpcp, selected, adj, params, semigreedy_vertex_selector);
    if (!success) {
      iterFile << ","
               << (col.get_n_colors() == 0
                       ? -1
                       : static_cast<int>(col.get_n_colors()));
      continue;
    } else {
      // Second step
      Col newCol;
      dpcp_dsatur_heur(dpcp, selected, adj, newCol);
      if (col.get_n_colors() == 0 ||
          newCol.get_n_colors() < col.get_n_colors()) {
        col = newCol;
        stats.bestTime =
            std::chrono::duration<double>(ClockType::now() - start).count();
        stats.bestIter = i;
      }
      iterFile << "," << col.get_n_colors();
    }
  }

  if (col.get_n_colors() == 0) {
    stats.state = UNKNOWN;
  } else {
    assert(col.check_coloring(dpcp));
    stats.state = FEASIBLE;
    stats.value = static_cast<double>(col.get_n_colors());
  }

  TimePoint end = ClockType::now();
  stats.totalTime = std::chrono::duration<double>(end - start).count();

  return stats;
}

/******************* */
/* 1-STEP HEURISTICS */
/******************* */

// Struct with the vertex information needed by the 1-step heuristic for
// DPCP
struct Heur1SVertexInfo {
  bool removed;               // whether the vertex has been removed
  Color color;                // color of the vertex (-1: uncolored)
  std::set<Color> adjColors;  // set of adjacent colors
  std::set<Vertex>
      uncolNeighbors;  // set of uncolored neighbors (outside of P[pi] U Q[qj])
  std::map<size_t, size_t> degree;  // Q-degree
  size_t pi;
  size_t qj;

  // Constructor
  Heur1SVertexInfo(const DPCPInst& dpcp, const Vertex& u)
      : removed(false),
        color(-1),
        adjColors(),
        uncolNeighbors(),
        pi(dpcp.get_P_part(u)),
        qj(dpcp.get_Q_part(u)) {
    for (Vertex v :
         boost::make_iterator_range(adjacent_vertices(u, dpcp.get_graph()))) {
      // Get components (pi,qj) of v
      size_t pi2 = dpcp.get_P_part(v);
      size_t qj2 = dpcp.get_Q_part(v);
      if (pi != pi2 && qj != qj2) {
        uncolNeighbors.insert(v);
        degree[qj2]++;
      }
    }
  }

  void print_info() {
    std::cout << "(pi,qj): (" << pi << "," << qj << "), color: " << color
              << ", adjColors: [";
    for (Color i : adjColors) std::cout << i << ",";
    std::cout << "], degree: " << uncolNeighbors.size()
              << ", removed: " << removed << std::endl;
  }
};

using InfoMap = std::map<Vertex, Heur1SVertexInfo>;
using PSizeMap = std::map<size_t, size_t>;

// Get the set of vertices that are invalidated when vertex v = (pi,qj) is
// colored with i, i.e. the returned set contains:
// (i) neighbors of v in Q[qj]
// (ii) vertices of Q[qj] that have i as an adjacent color
// (iii) neighbor (pi',qj') of v with pi != pi' and qj != qj' such
//       that some vertex of Q[qj'] has color i
// Raise an exception if v invalidates all remaining vertices of some P[pi']
VertexSet get_invalidated_vertices(const DPCPInst& dpcp, Col& col,
                                   const InfoMap& info, const PSizeMap& nP,
                                   Vertex v, Color i) {
  VertexSet inv;
  std::map<size_t, std::set<Vertex>> invMap;
  for (Vertex u : dpcp.get_Q()[dpcp.get_Q_part(v)])
    if (!info.at(u).removed && info.at(u).color == -1 &&
        ((!edge(u, v, dpcp.get_graph()).second &&
          info.at(u).adjColors.contains(i)) ||
         edge(u, v, dpcp.get_graph()).second)) {
      inv.insert(u);
      invMap[dpcp.get_P_part(u)].insert(u);
    }
  for (Vertex u : info.at(v).uncolNeighbors)
    if (col.is_colored_Q(dpcp.get_Q_part(u)) &&
        col.get_color_Q(dpcp.get_Q_part(u)) == i) {
      inv.insert(u);
      invMap[dpcp.get_P_part(u)].insert(u);
    }
  for (auto& [pi, vs] : invMap)
    if (vs.size() == nP.at(pi)) throw std::exception();
  return inv;
}

// Select the next vertex to color
// First, for each i in I, find the vertex with the minimum degree of
// saturation, where ties are broken by the minimum Q-degree.
// Seconds, among those vertices, choose the one with the maximum degree of
// saturation, where ties are broken by the minimum cardinality of P[pi]
Vertex greedy_vertex_selector_1S(const DPCPInst& dpcp, Col& col,
                                 const InfoMap& info, PSizeMap& nP,
                                 std::set<size_t>& unprocessedP) {
  // Choose a vertex with the maximum degree of saturation,
  // among the vertices with the minimum degree of saturation in each P[pi]
  // Ties are broken by the minimum cardinality of P_{i_v}
  // Further ties are broken by index
  Vertex maxv = NULL;
  size_t maxdsat = 0;
  size_t minPCardinality = std::numeric_limits<size_t>::max();

  for (auto pi : unprocessedP) {
    // Choose a vertex v in P[pi] with the minimum degree of saturation
    // Ties are broken by the minimum Q-degree, i.e., the number of Q_j
    // (other than j_v) that contains at least one uncolored neighbor of v
    Vertex minv = NULL;
    size_t mindsat = std::numeric_limits<size_t>::max();
    size_t minQdeg = std::numeric_limits<size_t>::max();
    for (Vertex v : dpcp.get_P()[pi]) {
      if (info.at(v).removed) continue;
      size_t dsat = info.at(v).adjColors.size();
      size_t Qdeg = info.at(v).degree.size();
      // Update the minimum vertex if needed
      if (minv == NULL || dsat < mindsat ||
          (dsat == mindsat && Qdeg < minQdeg)) {
        minv = v;
        mindsat = dsat;
        minQdeg = Qdeg;
      }
      // Early stop
      if (maxv != NULL && dsat < maxdsat) break;
    }
    // Update the maximum vertex if needed
    if (maxv == NULL || mindsat > maxdsat ||
        (mindsat == maxdsat && nP[pi] < minPCardinality)) {
      maxv = minv;
      maxdsat = mindsat;
      minPCardinality = nP[pi];
    }
  }

  return maxv;
}

// Select the color for vertex v
// When the color of v is not forced, the color is selected based on the
// minimum number of invalidated vertices. Colors that invalidate all
// remaining vertices in some P[pi'] are ignored. A new color is introduced
// only
// when mandatory (i.e. when the neighbors of v already use all available
// colors or the available colors force a failure)
std::pair<Color, VertexSet> greedy_color_selector_1S(const DPCPInst& dpcp,
                                                     Col& col,
                                                     const InfoMap& info,
                                                     PSizeMap& nP,
                                                     const Vertex v) {
  Color bestc = -1;
  VertexSet bestinv;

  // Is the color of v forced? That is, when Q[qj_v] has already colored
  // vertices
  if (col.is_colored_Q(dpcp.get_Q_part(v))) {
    bestc = col.get_color_Q(dpcp.get_Q_part(v));
    try {
      bestinv = get_invalidated_vertices(dpcp, col, info, nP, v, bestc);
    } catch (...) {
      throw std::exception();
    }
    return {bestc, bestinv};
  }

  // Otherwise, all colors not used by the neighbors of v are candidates
  // In addition, we consider a new color that increases the chromatic number
  std::list<Color> colors;
  for (size_t c = 0; c < col.get_n_colors(); ++c)
    if (!info.at(v).adjColors.contains(c)) colors.push_back(c);
  colors.push_back(col.get_n_colors());

  // Find a candidate color that invalidates the fewest number of vertices,
  // and ignoring colors that invalidate all remaining vertices in some P[pi']
  // Ties are broken by the color index
  for (auto c : colors) {
    bool isNewColor = (c == static_cast<int>(col.get_n_colors()));
    if (isNewColor && bestc != -1)
      continue;  // No need to try a new color if some used color is candidate
    VertexSet inv;
    try {
      inv = get_invalidated_vertices(dpcp, col, info, nP, v, c);
    } catch (...) {
      continue;
    }
    if (bestc == -1 || inv.size() < bestinv.size()) {
      bestc = c;
      bestinv = inv;
    }
  }

  if (bestc == -1) throw std::exception();

  return {bestc, bestinv};
}

// Single-step greedy heuristic for DPCP
// In each iteration, a vertex is selected and colored, with some criterion.
// Before proceeding to the next iteration, invalidated vertices are removed
bool single_step(const DPCPInst& dpcp, Col& col, bool greedy) {
  // Map with the necessary information of each vertex
  std::map<Vertex, Heur1SVertexInfo> info;

  // Map with the size of each P[pi]
  std::map<size_t, size_t> nP;
  for (size_t pi = 0; pi < dpcp.get_nP(); ++pi)
    nP.emplace(pi, dpcp.get_P()[pi].size());

  // Fill the information of the vertices and add them to the candidate list
  for (Vertex u : boost::make_iterator_range(vertices(dpcp.get_graph())))
    info.emplace(u, Heur1SVertexInfo(dpcp, u));

  // Subset of P indices that has not been processed yet
  std::set<size_t> unprocessedP;
  for (size_t pi = 0; pi < dpcp.get_nP(); ++pi) unprocessedP.insert(pi);

  while (!unprocessedP.empty()) {
    Vertex u = greedy_vertex_selector_1S(dpcp, col, info, nP, unprocessedP);
    std::pair<Color, VertexSet> tuple;
    try {
      tuple = greedy_color_selector_1S(dpcp, col, info, nP, u);
    } catch (...) {
      return false;
    }

    Color i = tuple.first;
    VertexSet& inv = tuple.second;
    unprocessedP.erase(info.at(u).pi);

    // Color u with i
    info.at(u).color = i;
    col.set_color(dpcp, dpcp.get_current_id(u), i);
    // info.at(u).print_info();

    // Invalidate the remaining vertices of P[pi_u]
    for (Vertex v : dpcp.get_P()[dpcp.get_P_part(u)])
      if (v != u && !info.at(v).removed) inv.insert(v);

    // Add i as an adjacent color
    for (Vertex v : info.at(u).uncolNeighbors) info.at(v).adjColors.insert(i);

    // Remove vertices and update information
    for (Vertex v : inv) {
      info.at(v).removed = true;
      nP.at(info.at(v).pi)--;
      for (Vertex w : info.at(v).uncolNeighbors) {
        info.at(w).uncolNeighbors.erase(v);
        size_t qjv = info.at(v).qj;
        info.at(w).degree[qjv]--;
        if (info.at(w).degree[qjv] == 0) info.at(w).degree.erase(qjv);
      }
    }
  }
  return true;
}

// 3-arg overload: discards iterFile output
HeurStats dpcp_2_step_semigreedy_heur(const DPCPInst& dpcp, Col& col,
                                      const Params& params) {
  struct NullBuffer : std::streambuf {
    int overflow(int c) { return c; }
  } nullBuffer;
  std::ostream nullstream(&nullBuffer);
  return dpcp_2_step_semigreedy_heur(dpcp, col, params, nullstream);
}

// One-step heuristic for DPCP
HeurStats dpcp_1_step_greedy_heur(const DPCPInst& dpcp, Col& col) {
  TimePoint start = ClockType::now();
  HeurStats stats;

  bool success = single_step(dpcp, col, true);

  if (!success) {
    TimePoint end = ClockType::now();
    stats.state = UNKNOWN;
    stats.totalTime = std::chrono::duration<double>(end - start).count();
    stats.totalIters = 1;
    return stats;
  }

  assert(col.check_coloring(dpcp));

  TimePoint end = ClockType::now();
  stats.state = FEASIBLE;
  stats.totalTime = std::chrono::duration<double>(end - start).count();
  stats.value = static_cast<double>(col.get_n_colors());
  stats.totalIters = 1;
  stats.bestTime = stats.totalTime;
  stats.bestIter = 0;

  return stats;
}
