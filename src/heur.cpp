#include "heur.hpp"
#include "graph.hpp"
#include "random.hpp"
extern "C" {
#include "color.h"
}

#include <chrono>
#include <limits>
#include <queue>
#include <random>
#include <ranges>
#include <stdexcept>

#define ALPHA_B 0.1 // Parameter for the semigreedy selection

using ClockType = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::_V2::system_clock::time_point;

/******************* */
/* 2-STEP HEURISTICS */
/******************* */

// Two-step greedy heuristic criteria

// Criteria 0: DEG-REAL
// Degree of v in G, ignoring neighbors in P_{i(v)} and Q_{j(v)}
size_t get_real_degree(const GraphEnv &genv,
                       const std::map<Vertex, bool> &removed,
                       const VertexVector &selected,
                       const std::map<TypeB, std::set<TypeB>> &adj,
                       const Vertex &v) {
  size_t degree = 0;
  for (Vertex u :
       boost::make_iterator_range(adjacent_vertices(v, genv.graph))) {
    if (removed.at(u))
      continue;
    if (genv.graph[u].first != genv.graph[v].first &&
        genv.graph[u].second != genv.graph[v].second)
      degree++;
  }
  return degree;
}

// Criteria 1: DEG-Q (used for breaking ties)
// Degree of v in G[Q_{j_v}]
size_t get_Q_degree(const GraphEnv &genv, const std::map<Vertex, bool> &removed,
                    const VertexVector &selected,
                    const std::map<TypeB, std::set<TypeB>> &adj,
                    const Vertex &v) {
  size_t neighbors = 0;
  for (Vertex u :
       boost::make_iterator_range(adjacent_vertices(v, genv.graph))) {
    if (removed.at(u))
      continue;
    if (genv.graph[u].second == genv.graph[v].second)
      neighbors++;
  }
  return neighbors;
}

// Criteria 2: DEG-COLLAPSED
// Degree of v in the collapsed graph according to Q,
// ignoring neighbors in P_{i(v)} and Q_{j(v)}
size_t get_collapsed_degree(const GraphEnv &genv,
                            const std::map<Vertex, bool> &removed,
                            const VertexVector &selected,
                            const std::map<TypeB, std::set<TypeB>> &adj,
                            const Vertex &v) {
  std::set<TypeB> Qdegree;
  for (Vertex u :
       boost::make_iterator_range(adjacent_vertices(v, genv.graph))) {
    if (removed.at(u))
      continue;
    if (genv.graph[u].first != genv.graph[v].first &&
        genv.graph[u].second != genv.graph[v].second)
      Qdegree.insert(genv.graph[u].second);
  }
  return Qdegree.size();
}

// Criteria 3: EDG
// No. of edges added to \tilde{G}[W]
size_t get_n_new_edge(const GraphEnv &genv,
                      const std::map<Vertex, bool> &removed,
                      const VertexVector &selected,
                      const std::map<TypeB, std::set<TypeB>> &adj,
                      const Vertex &v) {
  std::set<TypeB> adjBs;
  TypeB bv = genv.graph[v].second;
  for (Vertex u : selected) {
    if (!edge(v, u, genv.graph).second)
      continue;
    TypeB bu = genv.graph[u].second;
    if (bv == bu)
      continue;
    if (bv < bu && (!adj.contains(bv) || !adj.at(bv).contains(bu)))
      adjBs.insert(bu);
    else if (bv > bu && (!adj.contains(bu) || !adj.at(bu).contains(bv)))
      adjBs.insert(bu);
  }
  return adjBs.size();
}

size_t evaluate_vertex(const GraphEnv &genv,
                       const std::map<Vertex, bool> &removed,
                       const VertexVector &selected,
                       const std::map<TypeB, std::set<TypeB>> &adj,
                       size_t variant, const Vertex &v) {
  if (variant == 0)
    return get_real_degree(genv, removed, selected, adj, v);
  else if (variant == 1)
    return get_Q_degree(genv, removed, selected, adj, v);
  else if (variant == 2)
    return get_collapsed_degree(genv, removed, selected, adj, v);
  else
    return get_n_new_edge(genv, removed, selected, adj, v);
}

// Function to select the next vertex v following a greedy strategy
Vertex greedy_vertex_selector(const GraphEnv &genv,
                              const VertexVector &candidates,
                              const std::map<Vertex, bool> &removed,
                              const VertexVector &selected,
                              std::map<TypeB, std::set<TypeB>> &adj,
                              size_t variant) {
  Vertex minVertex = NULL;
  size_t minVal = std::numeric_limits<size_t>::max();
  size_t minTie = std::numeric_limits<size_t>::max();
  for (Vertex v : candidates) {
    if (removed.at(v))
      continue;
    size_t val = evaluate_vertex(genv, removed, selected, adj, variant, v);
    if (val < minVal) {
      minVal = val;
      minTie = evaluate_vertex(genv, removed, selected, adj, 1, v);
      minVertex = v;
    } else if (val == minVal) {
      size_t tie = evaluate_vertex(genv, removed, selected, adj, 1, v);
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
Vertex semigreedy_vertex_selector(const GraphEnv &genv,
                                  const VertexVector &candidates,
                                  const std::map<Vertex, bool> &removed,
                                  const VertexVector &selected,
                                  std::map<TypeB, std::set<TypeB>> &adj,
                                  size_t variant) {

  // First, find the lowest and highest degree among the candidates
  size_t minVal = std::numeric_limits<size_t>::max();
  size_t maxVal = 0;
  std::map<Vertex, size_t> valMap;
  for (Vertex v : candidates) {
    if (removed.at(v))
      continue;
    size_t val = evaluate_vertex(genv, removed, selected, adj, variant, v);
    valMap[v] = val;
    if (val < minVal)
      minVal = val;
    if (val > maxVal)
      maxVal = val;
  }
  double cutoff = minVal + ALPHA_B * (maxVal - minVal);

  // Build the RCL
  VertexVector rcl;
  for (auto [v, d] : valMap)
    if (d <= cutoff)
      rcl.push_back(v);

  if (rcl.empty())
    return NULL;

  // Select randomly a vertex from the RCL
  auto r = rand_int(rng) % rcl.size();

  return rcl[r];
}

using Heur2SVertexSelector = Vertex (*)(const GraphEnv &, const VertexVector &,
                                        const std::map<Vertex, bool> &,
                                        const VertexVector &,
                                        std::map<TypeB, std::set<TypeB>> &,
                                        size_t);

// Update the information after removing a vertex v
void update_info(const Graph &graph, std::map<Vertex, bool> &removed, Vertex v,
                 std::map<TypeA, size_t> &nVa) {
  removed.at(v) = true;
  TypeA av = graph[v].first;
  nVa.at(av)--;
}

// First step of the 2-step heuristic for DPCP
// Select a set of vertices of size |A|, one from each Va, such that the
// selected vertices of each Vb are a stable set
bool first_step(const GraphEnv &genv, VertexVector &selected,
                std::map<TypeB, std::set<TypeB>> &adj, size_t variant,
                Heur2SVertexSelector vertexSelector) {

  // Map with the removed vertices
  std::map<Vertex, bool> removed;
  for (Vertex u : boost::make_iterator_range(vertices(genv.graph)))
    removed.emplace(u, false);

  // Fill the map with the size of each Va
  std::map<TypeA, size_t> nVa;
  for (auto &[a, Va] : genv.Va)
    nVa.emplace(a, Va.size());

  // Subset of A that has not been processed yet
  std::set<TypeA> unprocessedA(genv.idA2TyA.begin(), genv.idA2TyA.end());

  // Comparison function for selecting the next a to process
  std::function<bool(const TypeA, const TypeA)> compareFunc;
  if (vertexSelector == greedy_vertex_selector) {
    // For greedy strategy, choose a such that Va has currently the minimum
    // size, break ties by the index of a
    compareFunc = [&nVa](const TypeA a1, const TypeA a2) {
      return (nVa.at(a1) < nVa.at(a2) || (nVa.at(a1) == nVa.at(a2) && a1 < a2));
    };
  } else {
    // Otherwise, choose a such that Va has currently the minimum size,
    // but break ties randomly
    compareFunc = [&nVa](const TypeA a1, const TypeA a2) {
      return (nVa.at(a1) < nVa.at(a2) ||
              (nVa.at(a1) == nVa.at(a2) && (rand_int(rng) % 2) == 0));
    };
  }

  while (selected.size() < genv.nA) {

    // Choose an unprocessed a
    TypeA a = *std::min_element(unprocessedA.begin(), unprocessedA.end(),
                                compareFunc);
    unprocessedA.erase(a);

    // Choose a vertex, with some criterion
    Vertex v =
        vertexSelector(genv, genv.Va.at(a), removed, selected, adj, variant);
    // info.at(v).print_info();
    if (v == NULL)
      return false; // No vertex can be selected from Va

    // Get label b of v
    TypeB b = genv.graph[v].second;

    // Update adjacent list and selected
    for (Vertex u : selected) {
      if (!edge(u, v, genv.graph).second)
        continue;
      TypeB b2 = genv.graph[u].second;
      if (b < b2)
        adj[b].insert(b2);
      else if (b > b2)
        adj[b2].insert(b);
    }
    selected.push_back(v);

    // Remove the vertices invalidated by choosing v, i.e.
    // vertices of Va \ v and neighbors of v that also belong to Vb
    for (Vertex u : genv.Va.at(a))
      if (!removed.at(u) && u != v)
        update_info(genv.graph, removed, u, nVa);
    for (Vertex u :
         boost::make_iterator_range(adjacent_vertices(v, genv.graph)))
      if (!removed.at(u) && genv.graph[u].second == b)
        update_info(genv.graph, removed, u, nVa);
  }

  return selected.size() == genv.nA;
}

// Second step of the 2-step heuristic for DPCP
// Color the selected vertices using DSATUR on the graph induced by them,
// collapsing vertices of each Vb
void dpcp_dsatur_heur(const GraphEnv &genv, VertexVector &selected,
                      std::map<TypeB, std::set<TypeB>> &adj, Col &col) {
  std::map<TypeB, size_t> bs; // Map from b to its new index (in the subgraph)
  std::vector<TypeB> invbs;   // Map from new index (in the subgraph) to b
  std::map<TypeB, VertexVector>
      repr;       // Map from b to the vector of vertices that b represents
  int ecount = 0; // Number of edges in the subgraph
  int elist[2 * num_edges(genv.graph)]; // Edgle list of the subgraph

  // First, initialize the bs, invbs and repr structures
  // This will be used for building the vertex set of the subgraph to color
  for (Vertex v : selected) {
    TypeB b1 = genv.graph[v].second;
    if (!bs.contains(b1)) {
      bs[b1] = bs.size();
      invbs.push_back(b1);
      repr[b1] = std::vector<Vertex>();
    }
    repr[b1].push_back(v);
  }

  // Seconds, build the edge list of the subgraph to color
  for (auto &[b1, adjB] : adj)
    for (TypeB b2 : adjB) {
      elist[2 * ecount] = bs[b1];
      elist[2 * ecount + 1] = bs[b2];
      ecount++;
    }

  // Third, color the subgraph using DSATUR
  int ncolors = 0;
  COLORset *colorclasses = NULL;
  COLORdsatur(bs.size(), ecount, elist, &ncolors, &colorclasses);

  // Recover coloring
  col.reset_coloring();
  for (int k = 0; k < ncolors; ++k)
    for (int j = 0; j < colorclasses[k].count; ++j) {
      TypeB b = invbs[colorclasses[k].members[j]];
      for (Vertex v : repr[b])
        col.set_color(genv.graph, v, k);
    }
}

// General two-step greedy heuristic for DPCP
Stats dpcp_2_step_greedy_heur(const GraphEnv &genv, Col &col, size_t variant) {

  TimePoint start = ClockType::now();
  Stats stats;

  // First step
  VertexVector selected;                // Vector of selected vertices
  std::map<TypeB, std::set<TypeB>> adj; // Adjacent list of the subgraph
                                        // induced by the selected vertices
  bool success =
      first_step(genv, selected, adj, variant, greedy_vertex_selector);
  if (!success) {
    stats.state = INFEASIBLE;
  } else {
    // Second step
    dpcp_dsatur_heur(genv, selected, adj, col);
    assert(col.check_coloring(genv.graph));
    stats.state = FEASIBLE;
    stats.ub = static_cast<double>(col.get_n_colors());
    stats.bestTime =
        std::chrono::duration<double>(ClockType::now() - start).count();
    stats.bestIter = 0;
  }

  TimePoint end = ClockType::now();
  stats.time = std::chrono::duration<double>(end - start).count();

  return stats;
}

// General two-step semigreedy heuristic for DPCP
Stats dpcp_2_step_semigreedy_heur(const GraphEnv &genv, Col &col, size_t nIters,
                                  size_t variant) {

  TimePoint start = ClockType::now();
  Stats stats;

  for (size_t i = 0; i < nIters; ++i) {

    // First step
    VertexVector selected;                // Vector of selected vertices
    std::map<TypeB, std::set<TypeB>> adj; // Adjacent list of the subgraph
                                          // induced by the selected vertices
    bool success =
        first_step(genv, selected, adj, variant, semigreedy_vertex_selector);
    if (!success) {
      continue;
    } else {
      // Second step
      Col newCol;
      dpcp_dsatur_heur(genv, selected, adj, newCol);
      if (col.get_n_colors() == 0 ||
          newCol.get_n_colors() < col.get_n_colors()) {
        col = newCol;
        stats.bestTime =
            std::chrono::duration<double>(ClockType::now() - start).count();
        stats.bestIter = i;
      }
    }
  }

  if (col.get_n_colors() == 0) {
    stats.state = INFEASIBLE;
  } else {
    assert(col.check_coloring(genv.graph));
    stats.state = FEASIBLE;
    stats.ub = static_cast<double>(col.get_n_colors());
  }

  TimePoint end = ClockType::now();
  stats.time = std::chrono::duration<double>(end - start).count();

  return stats;
}

/******************* */
/* 1-STEP HEURISTICS */
/******************* */

// Struct with the vertex information needed by the 1-step heuristic for
// DPCP
struct Heur1SVertexInfo {
  bool removed;              // whether the vertex has been removed
  Color color;               // color of the vertex (-1: uncolored)
  std::set<Color> adjColors; // set of adjacent colors
  std::set<Vertex>
      uncolNeighbors; // set of uncolored neighbors (outside of Va U Vb)
  std::map<TypeB, size_t> degree; // Q-degree
  TypeA a;
  TypeB b;

  // Constructor
  Heur1SVertexInfo(const Graph &graph, const Vertex &u)
      : removed(false), color(-1), adjColors(), uncolNeighbors(),
        a(graph[u].first), b(graph[u].second) {
    for (Vertex v : boost::make_iterator_range(adjacent_vertices(u, graph))) {
      // Get components (a,b) of v
      TypeA av = graph[v].first;
      TypeB bv = graph[v].second;
      if (a != av && b != bv) {
        uncolNeighbors.insert(v);
        degree[bv]++;
      }
    }
  }

  void print_info() {
    std::cout << "(a,b): (" << a << "," << b << "), color: " << color
              << ", adjColors: [";
    for (Color i : adjColors)
      std::cout << i << ",";
    std::cout << "], degree: " << uncolNeighbors.size()
              << ", removed: " << removed << std::endl;
  }
};

using InfoMap = std::map<Vertex, Heur1SVertexInfo>;
using VaSizeMap = std::map<TypeA, size_t>;

// Get the set of vertices that are invalidated when vertex v = (a,b) is
// colored with i, i.e. the returned set contains:
// (i) neighbors of v in V^b
// (ii) vertices of V^b that has i as an adjacent color
// (iii) neighbor (a',b') of v with a != a' and b != b' such
//       that some vertex of Vb' has color i
// Raise an exception if v invalidates all the remaining vertices of some V_a'
VertexSet get_invalidated_vertices(const GraphEnv &genv, Col &col,
                                   const InfoMap &info, const VaSizeMap &nVa,
                                   Vertex v, Color i) {
  VertexSet inv;
  std::map<TypeA, std::set<Vertex>> invMap;
  for (Vertex u : genv.Vb.at(genv.graph[v].second))
    if (!info.at(u).removed && info.at(u).color == -1 &&
        ((!edge(u, v, genv.graph).second && info.at(u).adjColors.contains(i)) ||
         edge(u, v, genv.graph).second)) {
      inv.insert(u);
      invMap[genv.graph[u].first].insert(u);
    }
  for (Vertex u : info.at(v).uncolNeighbors)
    if (col.is_colored_B(genv.graph[u].second) &&
        col.get_color_B(genv.graph[u].second) == i) {
      inv.insert(u);
      invMap[genv.graph[u].first].insert(u);
    }
  for (auto &[a, vs] : invMap)
    if (vs.size() == nVa.at(a))
      throw std::exception();
  return inv;
}

// Select the next vertex to color
// First, for each i in I, find the vertex with the minimum degree of
// saturation, where ties are broken by the minimum Q-degree.
// Seconds, among those vertices, choose the one with the maximum degree of
// saturation, where ties are broken by the minimum cardinality of P_i
Vertex greedy_vertex_selector_1S(const GraphEnv &genv, Col &col,
                                 const InfoMap &info, VaSizeMap &nPi, ASet &I) {

  // Choose a vertex with the maximum degree of saturation,
  // among the vertices with the minimum degree of saturation in each P_i
  // Ties are broken by the minimum cardinality of P_{i_v}
  // Further ties are broken by index
  Vertex maxv = NULL;
  size_t maxdsat = 0;
  size_t minnPi = std::numeric_limits<size_t>::max();

  for (auto i : I) {

    // Choose a vertex v in P_i with the minimum degree of saturation
    // Ties are broken by the minimum Q-degree, i.e., the number of Q_j
    // (other than j_v) that contains at least one uncolored neighbor of v
    Vertex minv = NULL;
    size_t mindsat = std::numeric_limits<size_t>::max();
    size_t minQdeg = std::numeric_limits<size_t>::max();
    for (Vertex v : genv.Va.at(i)) {
      if (info.at(v).removed)
        continue;
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
      if (maxv != NULL && dsat < maxdsat)
        break;
    }
    // Update the maximum vertex if needed
    if (maxv == NULL || mindsat > maxdsat ||
        (mindsat == maxdsat && nPi[i] < minnPi)) {
      maxv = minv;
      maxdsat = mindsat;
      minnPi = nPi[i];
    }
  }

  return maxv;
}

// Select the color for vertex v
// When the color of v is not forced, the color is selected based on the minimum
// number of invalidated vertices. Colors that invalidate all remaining vertices
// in some P_i' are ignored.
// A new color is introduced only when mandatory (i.e. when the neighbors of v
// already use all available colors or the available colors force a failure)
std::pair<Color, VertexSet>
greedy_color_selector_1S(const GraphEnv &genv, Col &col, const InfoMap &info,
                         VaSizeMap &nPi, const Vertex v) {

  Color bestc = -1;
  VertexSet bestinv;

  // Is the color of v forced? That is, when Q_{j_v} has already colored
  // vertices
  if (col.is_colored_B(genv.graph[v].second)) {
    bestc = col.get_color_B(genv.graph[v].second);
    try {
      bestinv = get_invalidated_vertices(genv, col, info, nPi, v, bestc);
    } catch (...) {
      throw std::exception();
    }
    return {bestc, bestinv};
  }

  // Otherwise, all colors not used by the neighbors of v are candidates
  // In addition, we consider a new color that increases the chromatic number
  std::list<Color> colors;
  for (size_t c = 0; c < col.get_n_colors(); ++c)
    if (!info.at(v).adjColors.contains(c))
      colors.push_back(c);
  colors.push_back(col.get_n_colors());

  // Find a candidate color that invalidates the fewest number of vertices,
  // and ignoring colors that invalidate all remaining vertices in some P_i'
  // Ties are broken by the color index
  for (auto c : colors) {
    bool isNewColor = (c == static_cast<int>(col.get_n_colors()));
    if (isNewColor && bestc != -1)
      continue; // No need to try a new color if some used color is candidate
    VertexSet inv;
    try {
      inv = get_invalidated_vertices(genv, col, info, nPi, v, c);
    } catch (...) {
      continue;
    }
    if (bestc == -1 || inv.size() < bestinv.size()) {
      bestc = c;
      bestinv = inv;
    }
  }

  if (bestc == -1)
    throw std::exception();

  return {bestc, bestinv};
}

// Single-step greedy heuristic for DPCP
// In each iteration, a vertex is selected and colored, with some criterion.
// Before proceeding to the next iteration, invalidated vertices are removed
bool single_step(const GraphEnv &genv, Col &col, bool greedy) {

  // Map with the necessary information of each vertex
  std::map<Vertex, Heur1SVertexInfo> info;

  // Map with the size of each Va
  std::map<TypeA, size_t> nVa;
  for (auto &[a, Va] : genv.Va)
    nVa.emplace(a, Va.size());

  // Fill the information of the vertices and add them to the candidate list
  for (Vertex u : boost::make_iterator_range(vertices(genv.graph)))
    info.emplace(u, Heur1SVertexInfo(genv.graph, u));

  // Subset of A that has not been processed yet
  ASet A(genv.idA2TyA.begin(), genv.idA2TyA.end());

  while (!A.empty()) {

    Vertex u = greedy_vertex_selector_1S(genv, col, info, nVa, A);
    std::pair<Color, VertexSet> tuple;
    try {
      tuple = greedy_color_selector_1S(genv, col, info, nVa, u);
    } catch (...) {
      return false;
    }

    Color i = tuple.first;
    VertexSet &inv = tuple.second;
    A.erase(info.at(u).a);

    // Color u with i
    info.at(u).color = i;
    col.set_color(genv.graph, u, i);
    // info.at(u).print_info();

    // Invalidate the remaining vertices of Va
    for (Vertex v : genv.Va.at(genv.graph[u].first))
      if (v != u && !info.at(v).removed)
        inv.insert(v);

    // Add i as an adjacent color
    for (Vertex v : info.at(u).uncolNeighbors)
      info.at(v).adjColors.insert(i);

    // Remove vertices and update information
    for (Vertex v : inv) {
      info.at(v).removed = true;
      nVa.at(info.at(v).a)--;
      for (Vertex w : info.at(v).uncolNeighbors) {
        info.at(w).uncolNeighbors.erase(v);
        TypeB jv = info.at(v).b;
        info.at(w).degree[jv]--;
        if (info.at(w).degree[jv] == 0)
          info.at(w).degree.erase(jv);
      }
    }
  }
  return true;
}

// One-step heuristic for DPCP
Stats dpcp_1_step_greedy_heur(const GraphEnv &genv, Col &col) {
  TimePoint start = ClockType::now();
  Stats stats;

  bool success = single_step(genv, col, true);

  if (!success) {
    TimePoint end = ClockType::now();
    stats.state = INFEASIBLE;
    stats.time = std::chrono::duration<double>(end - start).count();
    return stats;
  }

  assert(col.check_coloring(genv.graph));

  TimePoint end = ClockType::now();
  stats.state = FEASIBLE;
  stats.time = std::chrono::duration<double>(end - start).count();
  stats.ub = static_cast<double>(col.get_n_colors());
  stats.bestTime = stats.time;
  stats.bestIter = 0;

  return stats;
}
