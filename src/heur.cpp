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

// Get the current degree of v = (a,b) in Vb, ignoring removed vertices
size_t get_current_degree_in_Vb(const GraphEnv &genv,
                                const std::map<Vertex, bool> &removed,
                                const VertexVector &selected,
                                const std::map<TypeB, std::set<TypeB>> &adj,
                                const Vertex &v) {
  size_t degree = 0;
  for (Vertex u : boost::make_iterator_range(adjacent_vertices(v, genv.graph)))
    if (!removed.at(u) && genv.graph[u].second == genv.graph[v].second)
      degree++;
  return degree;
}

// Get the current number of adjacencies of v = (a,b) in B, ignoring removed
// vertices, i.e the size of {b' : (a',b') is not removed, (a,b) is
// adjacent to (a',b')}
size_t get_current_degree_in_B(const GraphEnv &genv,
                               const std::map<Vertex, bool> &removed,
                               const VertexVector &selected,
                               const std::map<TypeB, std::set<TypeB>> &adj,
                               const Vertex &v) {
  std::set<TypeB> adjBs;
  for (Vertex u : boost::make_iterator_range(adjacent_vertices(v, genv.graph)))
    if (!removed.at(u) && genv.graph[u].second != genv.graph[v].second)
      adjBs.insert(genv.graph[u].second);
  return adjBs.size();
}

// Get the current number of adjacencies of v = (a,b) in B, restricted to the
// selected vertices, i.e the size of {b' : (a',b') is selected, (a,b) is
// adjacent to (a',b')}
size_t get_current_degree_in_selected_B(
    const GraphEnv &genv, const std::map<Vertex, bool> &removed,
    const VertexVector &selected, const std::map<TypeB, std::set<TypeB>> &adj,
    const Vertex &v) {
  std::set<TypeB> adjBs;
  for (Vertex u : selected)
    if (edge(v, u, genv.graph).second &&
        genv.graph[u].second != genv.graph[v].second)
      adjBs.insert(genv.graph[u].second);
  return adjBs.size();
}

// Get the number of new edges that would be added after selecting v = (a,b)
// A new edge is added when there is a selected vertex u = (a',b') such that
// u is adjacent to v and there is no edge bb' yet
size_t get_n_new_edges(const GraphEnv &genv,
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
    return get_current_degree_in_Vb(genv, removed, selected, adj, v);
  else if (variant == 1)
    return get_current_degree_in_B(genv, removed, selected, adj, v);
  else if (variant == 2)
    return get_current_degree_in_selected_B(genv, removed, selected, adj, v);
  else
    return get_n_new_edges(genv, removed, selected, adj, v);
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
  for (Vertex v : candidates) {
    if (removed.at(v))
      continue;
    size_t val = evaluate_vertex(genv, removed, selected, adj, variant, v);
    if (val < minVal) {
      minVal = val;
      minVertex = v;
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
void second_step(const GraphEnv &genv, VertexVector &selected,
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
    second_step(genv, selected, adj, col);
    assert(col.check_coloring(genv.graph));
    stats.state = FEASIBLE;
    stats.ub = static_cast<double>(col.get_n_colors());
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

  while (nIters-- > 0) {

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
      second_step(genv, selected, adj, newCol);
      if (col.get_n_colors() == 0 || newCol.get_n_colors() < col.get_n_colors())
        col = newCol;
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
      if (a != av && b != bv)
        uncolNeighbors.insert(v);
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

// Select a set of candidate vertices to color, with the following criterion:
// (i) Among the uncolored a, choose one with the lowest |V_a|
// (ii) Build a set of candidates, containing every vertex (a,b) of V_a
// such that: (1) some vertex in V^b is already colored or (2) (a,b) has
// the lowest degree of saturation in V_a
VertexVector vertices_selector(const GraphEnv &genv, Col &col,
                               const InfoMap &info, VaSizeMap &nVa, ASet &A,
                               bool greedy) {

  // Comparison function for elements of A
  std::function<bool(const TypeA, const TypeA)> compareFunc;
  if (greedy) {
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

  // Choose the next a
  TypeA a = *std::min_element(A.begin(), A.end(), compareFunc);
  A.erase(a);

  // Set of candidates
  VertexVector candidates1, candidates2;
  size_t minDSat = std::numeric_limits<size_t>::max();
  for (Vertex v : genv.Va.at(a)) {
    if (info.at(v).removed)
      continue;
    TypeB b = genv.graph[v].second;
    if (col.is_colored_B(b)) {
      candidates1.push_back(v);
      continue;
    }
    size_t dsat = info.at(v).adjColors.size();
    if (dsat < minDSat) {
      minDSat = dsat;
      candidates2.clear();
    }
    if (dsat == minDSat)
      candidates2.push_back(v);
  }
  for (Vertex v : candidates1)
    candidates2.push_back(v);

  return candidates2;
}

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

// Select a vertex from the candidates and a used free color, with the following
// criterion: the one that invalidates the lowest number of vertices.
// If there is none, do the same with a new color.
// Return the vertex, the color, and the set of invalidated vertices
std::tuple<Vertex, Color, VertexSet>
vertex_color_selector(const GraphEnv &genv, Col &col, const InfoMap &info,
                      const VaSizeMap &nVa, const VertexVector &vCandidates) {

  VertexSet minInv;
  size_t n_minInv = std::numeric_limits<size_t>::max();
  Vertex v = NULL;
  Color c = -1;

  // Iterate over the candidate vertices
  for (Vertex u : vCandidates) {

    TypeA a = genv.graph[u].first;
    TypeB b = genv.graph[u].second;

    // Find candidate colors for u
    std::vector<Color> colors;
    if (col.is_colored_B(b))
      colors.push_back(col.get_color_B(b));
    else {
      // A new color is needed? Ignore!
      if (info.at(u).adjColors.size() == col.get_n_colors())
        continue;
      // Otherwse, find the non-adjacent colors
      else {
        for (size_t i = 0; i < col.get_n_colors(); ++i)
          if (!info.at(u).adjColors.contains(i))
            colors.push_back(i);
      }
    }

    // Among the "safe" vertex-color pairs, i.e. that do not invalidate all the
    // remaining vertices of any V_a', select the one that minimice the number
    // of invalidated vertices
    for (Color k : colors) {
      VertexSet inv;
      try {
        inv = get_invalidated_vertices(genv, col, info, nVa, u, k);
      } catch (...) {
        continue;
      }
      if (inv.size() < n_minInv) {
        n_minInv = inv.size();
        minInv = inv;
        v = u;
        c = k;
      }
    }
  }

  // If any pair is candidate, find the new color that minimice the number
  // of invalidated vertices
  if (v == NULL) {

    // Iterate over the candidate vertices
    for (Vertex u : vCandidates) {

      // Ignore vertices with fixed color
      if (col.is_colored_B(genv.graph[u].second))
        continue;

      VertexSet inv;
      try {
        inv = get_invalidated_vertices(genv, col, info, nVa, u,
                                       col.get_n_colors());
      } catch (...) {
        continue;
      }
      if (inv.size() < n_minInv) {
        n_minInv = inv.size();
        minInv = inv;
        v = u;
        c = col.get_n_colors();
      }
    }
  }

  // If any pair is candidate yet, throw exception
  if (v == NULL)
    throw std::exception();

  return std::make_tuple(v, c, minInv);
}

// One-step heuristic for DPCP
// Based on a DSATUR implementation for GCP that runs in O(n^2)
Stats dpcp_1_step_greedy_heur(const GraphEnv &genv, Col &col) {
  TimePoint start = ClockType::now();

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

    auto vCandidates = vertices_selector(genv, col, info, nVa, A, true);
    std::tuple<Vertex, Color, VertexSet> tuple;
    try {
      tuple = vertex_color_selector(genv, col, info, nVa, vCandidates);
    } catch (...) {
      TimePoint end = ClockType::now();
      Stats stats;
      stats.state = INFEASIBLE;
      stats.time = std::chrono::duration<double>(end - start).count();
      return stats;
    }

    Vertex u = std::get<0>(tuple);
    Color i = std::get<1>(tuple);
    VertexSet inv = std::get<2>(tuple);

    // Get components (a,b) of u
    TypeA a = genv.graph[u].first;
    TypeB b = genv.graph[u].second;

    // Color u with i
    // std::cout << "Pintando: (" << genv.graph[u].first << ","
    //           << genv.graph[u].second << ") con " << i << std::endl;
    info.at(u).color = i;
    col.set_color(genv.graph, u, i);

    // Invalidate the remaining vertices of Va
    for (Vertex v : genv.Va.at(a))
      if (v != u && !info.at(v).removed)
        inv.insert(v);

    // Add i as an adjacent color
    for (Vertex v : info.at(u).uncolNeighbors)
      info.at(v).adjColors.insert(i);

    // Remove vertices and update information
    for (Vertex v : inv) {
      info.at(v).removed = true;
      TypeA av = info.at(v).a;
      nVa.at(av)--;
      for (Vertex w : info.at(v).uncolNeighbors)
        info.at(w).uncolNeighbors.erase(v);
    }
  }

  assert(col.check_coloring(genv.graph));

  TimePoint end = ClockType::now();

  Stats stats;
  stats.state = FEASIBLE;
  stats.time = std::chrono::duration<double>(end - start).count();
  stats.ub = static_cast<double>(col.get_n_colors());

  return stats;
}
