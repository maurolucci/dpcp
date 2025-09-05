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

#define ALPHA_B 0.1 // Parameter for the semigreedy selection

using ClockType = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::_V2::system_clock::time_point;

/******************* */
/* 2-STEP HEURISTICS */
/******************* */

using Heur2SDegreeFunc = size_t (*)(const GraphEnv &,
                                    const std::map<Vertex, bool> &,
                                    const VertexVector &,
                                    const std::map<TypeB, std::set<TypeB>> &,
                                    const Vertex &);

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

using Heur2SVertexSelector = Vertex (*)(const GraphEnv &, const VertexVector &,
                                        const std::map<Vertex, bool> &,
                                        const VertexVector &,
                                        std::map<TypeB, std::set<TypeB>> &,
                                        Heur2SDegreeFunc);

// Function to select the next vertex v = (a,b) such that v has the minimum
// degree in Vb following a greedy strategy
Vertex greedy_vertex_selector(const GraphEnv &genv,
                              const VertexVector &candidates,
                              const std::map<Vertex, bool> &removed,
                              const VertexVector &selected,
                              std::map<TypeB, std::set<TypeB>> &adj,
                              Heur2SDegreeFunc degreeFunc) {
  Vertex minVertex = NULL;
  size_t minDegree = std::numeric_limits<size_t>::max();
  for (Vertex v : candidates) {
    if (removed.at(v))
      continue;
    size_t currentDegree = degreeFunc(genv, removed, selected, adj, v);
    if (currentDegree < minDegree) {
      minDegree = currentDegree;
      minVertex = v;
    }
  }
  return minVertex;
}

// Function to select the next vertex v = (a,b) such that v has the minimum
// degree in Vb following a semigreedy strategy
// We use a restricted candidate list (RCL) of size |RCL| = max(1, alpha*|C|)
// where C is the set of candidates and alpha in [0,1] is a parameter
// We store in the RCL the best candidates and we select randomly one of them
Vertex semigreedy_vertex_selector(const GraphEnv &genv,
                                  const VertexVector &candidates,
                                  const std::map<Vertex, bool> &removed,
                                  const VertexVector &selected,
                                  std::map<TypeB, std::set<TypeB>> &adj,
                                  Heur2SDegreeFunc degreeFunc) {

  // First, find the lowest and highest degree among the candidates
  size_t minDegree = std::numeric_limits<size_t>::max();
  size_t maxDegree = 0;
  std::map<Vertex, size_t> degreeMap;
  for (Vertex v : candidates) {
    if (removed.at(v))
      continue;
    size_t currentDegree = degreeFunc(genv, removed, selected, adj, v);
    degreeMap[v] = currentDegree;
    if (currentDegree < minDegree)
      minDegree = currentDegree;
    if (currentDegree > maxDegree)
      maxDegree = currentDegree;
  }
  double cutoff = minDegree + ALPHA_B * (maxDegree - minDegree);

  // Build the RCL
  VertexVector rcl;
  for (auto [v, d] : degreeMap)
    if (d <= cutoff)
      rcl.push_back(v);

  if (rcl.empty())
    return NULL;

  // Select randomly a vertex from the RCL
  auto r = rand_int(rng) % rcl.size();

  return rcl[r];
}

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
                std::map<TypeB, std::set<TypeB>> &adj,
                Heur2SVertexSelector vertexSelector,
                Heur2SDegreeFunc degreeFunc) {

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
        vertexSelector(genv, genv.Va.at(a), removed, selected, adj, degreeFunc);
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
Stats dpcp_2_step_greedy_heur(const GraphEnv &genv, Col &col,
                              Heur2SDegreeFunc degreeFunc) {

  TimePoint start = ClockType::now();
  Stats stats;

  // First step
  VertexVector selected;                // Vector of selected vertices
  std::map<TypeB, std::set<TypeB>> adj; // Adjacent list of the subgraph
                                        // induced by the selected vertices
  bool success =
      first_step(genv, selected, adj, greedy_vertex_selector, degreeFunc);
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
                                  Heur2SDegreeFunc degreeFunc) {

  TimePoint start = ClockType::now();
  Stats stats;

  while (nIters-- > 0) {

    // First step
    VertexVector selected;                // Vector of selected vertices
    std::map<TypeB, std::set<TypeB>> adj; // Adjacent list of the subgraph
                                          // induced by the selected vertices
    bool success =
        first_step(genv, selected, adj, semigreedy_vertex_selector, degreeFunc);
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

Stats dpcp_2_step_greedy_heur_min_degree_in_Vb(const GraphEnv &genv, Col &col) {
  return dpcp_2_step_greedy_heur(genv, col, get_current_degree_in_Vb);
}
Stats dpcp_2_step_greedy_heur_min_degree_in_B(const GraphEnv &genv, Col &col) {
  return dpcp_2_step_greedy_heur(genv, col, get_current_degree_in_B);
}
Stats dpcp_2_step_greedy_heur_min_degree_in_selected_B(const GraphEnv &genv,
                                                       Col &col) {
  return dpcp_2_step_greedy_heur(genv, col, get_current_degree_in_selected_B);
}
Stats dpcp_2_step_greedy_heur_min_n_new_edges(const GraphEnv &genv, Col &col) {
  return dpcp_2_step_greedy_heur(genv, col, get_n_new_edges);
}

Stats dpcp_2_step_semigreedy_heur_min_degree_in_Vb(const GraphEnv &genv,
                                                   Col &col, size_t nIters) {
  return dpcp_2_step_semigreedy_heur(genv, col, nIters,
                                     get_current_degree_in_Vb);
}
Stats dpcp_2_step_semigreedy_heur_min_degree_in_B(const GraphEnv &genv,
                                                  Col &col, size_t nIters) {
  return dpcp_2_step_semigreedy_heur(genv, col, nIters,
                                     get_current_degree_in_B);
}
Stats dpcp_2_step_semigreedy_heur_min_degree_in_selected_B(const GraphEnv &genv,
                                                           Col &col,
                                                           size_t nIters) {
  return dpcp_2_step_semigreedy_heur(genv, col, nIters,
                                     get_current_degree_in_selected_B);
}
Stats dpcp_2_step_semigreedy_heur_min_n_new_edges(const GraphEnv &genv,
                                                  Col &col, size_t nIters) {
  return dpcp_2_step_semigreedy_heur(genv, col, nIters, get_n_new_edges);
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

// Given v = (a,bv) and u = (a,bu), v < u if (lexicographic order)
// (i) v has a lower saturation degree
// (ii) v has a lower uncolored degree
// (iii) bv has a lower index
bool vertex_comparison_VC1_1(const Heur1SVertexInfo &v,
                             const Heur1SVertexInfo &u) {
  assert(v.a == u.a);
  if (v.adjColors.size() < u.adjColors.size() ||
      (v.adjColors.size() == u.adjColors.size() &&
       v.uncolNeighbors.size() < u.uncolNeighbors.size()) ||
      (v.adjColors.size() == u.adjColors.size() &&
       v.uncolNeighbors.size() == u.uncolNeighbors.size() && v.b < u.b))
    return true;
  return false;
}

// Given v = (av,bv) and u = (au,bu) with av != au, v < u if (lexicographic
// order)
// (i) v has a greater saturation degree
// (ii) v has a greater uncolored degree
// (iii) av has a lower index
bool vertex_comparison_VC1_2(const Heur1SVertexInfo &v,
                             const Heur1SVertexInfo &u) {
  assert(v.a != u.a);
  if (v.adjColors.size() > u.adjColors.size() ||
      (v.adjColors.size() == u.adjColors.size() &&
       v.uncolNeighbors.size() > u.uncolNeighbors.size()) ||
      (v.adjColors.size() == u.adjColors.size() &&
       v.uncolNeighbors.size() == u.uncolNeighbors.size() && v.a < u.a))
    return true;
  return false;
}

// Function to select the next vertex to color using the vertex criterion
// VC1
Vertex vertex_selector_VC1(const std::set<Vertex> &candidates,
                           const std::map<Vertex, Heur1SVertexInfo> &info) {
  std::map<TypeA, Vertex> best;
  for (Vertex v : candidates) {
    TypeA a = info.at(v).a;
    if (!best.contains(a))
      best.insert(std::make_pair(a, v));
    else if (vertex_comparison_VC1_1(info.at(v), info.at(best.at(a))))
      best.at(a) = v;
  }
  return std::min_element(best.begin(), best.end(),
                          [&info](auto p1, auto p2) {
                            return vertex_comparison_VC1_2(info.at(p1.second),
                                                           info.at(p2.second));
                          })
      ->second;
}

// Check if painting a vertex v = (a,b) with color i is safe,
// i.e. it does not invalidate any vertex apart from Va and [Vb cap N(v)]
bool is_safe_color(const GraphEnv &genv, Col &col,
                   const std::map<Vertex, Heur1SVertexInfo> &info, Vertex v,
                   Color i) {
  for (Vertex u : genv.Vb.at(genv.graph[v].second))
    if (!info.at(u).removed && !edge(u, v, genv.graph).second &&
        info.at(u).adjColors.contains(i))
      return false;
  for (Vertex u : info.at(v).uncolNeighbors)
    if (col.is_colored_B(genv.graph[u].second) &&
        col.get_color_B(genv.graph[u].second) == i)
      return false;
  return true;
}

// Function to select a color for a vertex using the color criterion CC1
// Criterion CC1: assign the first free and *safe* color
Color color_selector_CC1(const GraphEnv &genv, Col &col,
                         const std::map<Vertex, Heur1SVertexInfo> &info,
                         const Vertex u) {
  // A new color is needed?
  if (info.at(u).adjColors.size() == col.get_n_colors())
    return col.get_n_colors();

  // Mark the colors used in the neighborhood of u
  std::vector<bool> used(std::min(genv.nA, genv.nB));
  for (Color i : info.at(u).adjColors)
    used[i] = true;

  // Find the first available and safe color for u
  Color i;
  for (i = 0; i < static_cast<int>(used.size()); i++)
    if (used[i] == false && is_safe_color(genv, col, info, u, i))
      break;

  return i;
}

// Get the number of vertices that are invalidated when vertex v = (a,b) is
// painted with color i
bool n_invalidated_vertices(const GraphEnv &genv, Col &col,
                            const std::map<Vertex, Heur1SVertexInfo> &info,
                            Vertex v, Color i) {
  size_t nInv = 0;
  for (Vertex u : genv.Vb.at(genv.graph[v].second))
    if (!info.at(u).removed && !edge(u, v, genv.graph).second &&
        info.at(u).adjColors.contains(i))
      nInv++;
  for (Vertex u : info.at(v).uncolNeighbors)
    if (col.is_colored_B(genv.graph[u].second) &&
        col.get_color_B(genv.graph[u].second) == i)
      nInv++;
  return nInv;
}

// Function to select a color for a vertex using the color criterion CC2
// Criterion CC2: between all free colors, assign the one that invalidates
// the lowest number of vertices
// TODO: ignore colors that force infeasibility, i.e, those colors that
// invalidate all the remaining vertices of some Va
Color color_selector_CC2(const GraphEnv &genv, Col &col,
                         const std::map<Vertex, Heur1SVertexInfo> &info,
                         const Vertex u) {
  // A new color is needed?
  if (info.at(u).adjColors.size() == col.get_n_colors())
    return col.get_n_colors();

  // Set of avilable colors (not used in the neighbors)
  std::set<Color> freeColors;
  for (size_t i = 0; i < col.get_n_colors(); ++i)
    if (!info.at(u).adjColors.contains(i))
      freeColors.insert(i);

  // Find the available color that invalidates the lowest number of vertices
  return *std::min_element(
      freeColors.begin(), freeColors.end(),
      [&genv, &col, &info, &u](const Color &i, const Color &j) {
        return n_invalidated_vertices(genv, col, info, u, i) <
               n_invalidated_vertices(genv, col, info, u, j);
      });
}

// One-step heuristic for DPCP
// Based on a DSATUR implementation for GCP that runs in O(n^2)
Stats dpcp_heur_1_step(const GraphEnv &genv, Col &col) {
  TimePoint start = ClockType::now();

  // Map with the necessary information of each vertex
  std::map<Vertex, Heur1SVertexInfo> info;

  // Map with the size of each Va
  std::map<TypeA, size_t> nVa;
  for (auto &[a, Va] : genv.Va)
    nVa.emplace(a, Va.size());

  // Set of candidates vertices
  std::set<Vertex> candidates;

  // Fill the information of the vertices and add them to the candidate list
  for (Vertex u : boost::make_iterator_range(vertices(genv.graph))) {
    info.emplace(u, Heur1SVertexInfo(genv.graph, u));
    candidates.insert(u);
  }

  while (!candidates.empty()) {

    // Choose the next vertex
    Vertex u = vertex_selector_VC1(candidates, info);
    info.at(u).print_info();

    // Get components (a,b) of u
    TypeA a = genv.graph[u].first;
    TypeB b = genv.graph[u].second;

    // Decide the color of u
    // Recall: all colored vertices in V^b must have the same color
    Color i = col.is_colored_B(b) ? col.get_color_B(b)
                                  : color_selector_CC2(genv, col, info, u);

    // Color u with i
    // std::cout << "Pintando: (" << genv.graph[u].first << ","
    //           << genv.graph[u].second << ") con " << i << std::endl;
    info.at(u).color = i;
    col.set_color(genv.graph, u, i);

    // Set of vertices that must be removed after coloring u
    std::set<Vertex> toRemove;

    // 1. Remove u and neighbors of u in Va
    // Recall: Va is a clique
    for (Vertex v : genv.Va.at(a))
      if (!info.at(v).removed)
        toRemove.insert(v);

    // 2. Remove neighbors of u in Vb
    // 3. Remove the vertices of V^b that has i as an adjacent color
    for (Vertex v : genv.Vb.at(b))
      if (!info.at(v).removed && info.at(v).color == -1 &&
          (edge(u, v, genv.graph).second || info.at(v).adjColors.contains(i)))
        toRemove.insert(v);

    // 4. Remove every neighbor (a',b') of u with a != a' and b != b' such
    // that some vertex of Vb' has color i Also, add i as an adjacent color
    for (Vertex v : info.at(u).uncolNeighbors) {
      // Get component b of v
      TypeB bv = genv.graph[v].second;
      if (col.is_colored_B(bv) && col.get_color_B(bv) == i)
        toRemove.insert(v);
      else
        info.at(v).adjColors.insert(i);
    }

    // Remove vertices and update information
    // Also, check for empty Va
    for (Vertex v : toRemove) {
      candidates.erase(v);
      info.at(v).removed = false;
      TypeA av = info.at(v).a;
      nVa.at(av)--;
      if (av != a && nVa.at(av) == 0) {
        // Cannot find solution
        TimePoint end = ClockType::now();
        Stats stats;
        stats.state = INFEASIBLE;
        stats.time = std::chrono::duration<double>(end - start).count();
        return stats;
      }
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

bool is_conflict(const GraphEnv &genv, Vertex v, Vertex u) {
  return edge(v, u, genv.graph).second;
}

void heur_solve_aux(const GraphEnv &genv, const std::vector<TypeA> &as,
                    Col &col) {
  VertexVector vertices;      // Vector of chosen vertices to be colored
  std::map<TypeB, size_t> bs; // Map from b to its new index (in the subgraph)
  std::vector<TypeB> invbs;   // Map from new index (in the subgraph) to b
  std::map<TypeB, VertexVector>
      repr; // Map from b to the vector of vertices that b represents
  std::vector<std::vector<TypeB>> adj;  // Adjacent lists of the subgraph
                                        // : new index -> std::vector<TypeB>
  int ecount = 0;                       // Number of edges in the subgraph
  int elist[2 * num_edges(genv.graph)]; // Edgle list of the subgraph

  col.reset_coloring();

  for (TypeA a : as) {
    // For each a in A, choose v in snd[i_a] such that:
    // v minimizes the conflicts with the already chosen representatives
    size_t minNNewConflicts = std::numeric_limits<size_t>::max();
    // size_t minNAllConflicts = std::numeric_limits<size_t>::max();
    std::set<TypeB> minNewConflicts;
    // std::set<TypeB> minAllConflicts;
    Vertex bestVertex;
    for (Vertex v : genv.snd[genv.tyA2idA.at(a)]) {
      std::set<TypeB> newConflicts;
      std::set<TypeB> allConflicts;
      TypeB bv = genv.graph[v].second;
      for (Vertex u : vertices) {
        TypeB bu = genv.graph[u].second;
        if (is_conflict(genv, v, u)) {
          if (!bs.contains(bv) ||
              find(adj[bs.at(bv)].begin(), adj[bs.at(bv)].end(), bu) ==
                  adj[bs[bv]].end())
            newConflicts.insert(bu);
          // allConflicts.insert(bu);
        }
      }
      if ((newConflicts.size() < minNNewConflicts) ||
          (newConflicts.size() == minNNewConflicts && rand_double(rng) < 0.5)) {
        minNNewConflicts = newConflicts.size();
        // minNAllConflicts = allConflicts.size();
        // minAllConflicts = allConflicts;
        minNewConflicts = newConflicts;
        bestVertex = v;
      }
    }

    vertices.push_back(bestVertex);
    TypeB b1 = genv.graph[bestVertex].second;
    if (!bs.contains(b1)) {
      bs[b1] = bs.size();
      invbs.push_back(b1);
      adj.push_back(std::vector<TypeB>());
      repr[b1] = std::vector<Vertex>();
    }
    repr[b1].push_back(bestVertex);
    for (auto b2 : minNewConflicts) {
      elist[2 * ecount] = bs[b1];
      elist[2 * ecount++ + 1] = bs[b2];
      adj[bs[b1]].push_back(b2);
      adj[bs[b2]].push_back(b1);
    }
  }

  int ncolors = 0;
  COLORset *colorclasses = NULL;
  COLORdsatur(bs.size(), ecount, elist, &ncolors, &colorclasses);

  // Recover coloring
  for (int k = 0; k < ncolors; ++k)
    for (int j = 0; j < colorclasses[k].count; ++j) {
      TypeB b = invbs[colorclasses[k].members[j]];
      for (Vertex v : repr[b])
        col.set_color(genv.graph, v, k);
    }
  // assert(col.check_coloring());
  return;
}

Stats heur_solve(const GraphEnv &genv, const std::vector<TypeA> &as, Col &col,
                 size_t repetitions, Pool &pool) {
  TimePoint start = ClockType::now();
  Col currentCol;
  std::vector<TypeA> ascopy(as);
  size_t minNCol = std::numeric_limits<size_t>::max();
  for (size_t iter = 0; iter < repetitions; iter++) {
    std::shuffle(ascopy.begin(), ascopy.end(), rng);
    heur_solve_aux(genv, ascopy, currentCol);
    if (currentCol.get_n_colors() < minNCol) {
      minNCol = currentCol.get_n_colors();
      col = currentCol;
    }
    // Add the stable to the pool
    for (size_t k = 0; k < currentCol.get_n_colors(); ++k)
      pool.push_back(currentCol.get_stable(genv.graph, k));
  }
  assert(col.check_coloring(genv.graph));
  TimePoint end = ClockType::now();

  Stats stats;
  stats.state = FEASIBLE;
  stats.time = std::chrono::duration<double>(end - start).count();
  stats.ub = static_cast<double>(col.get_n_colors());

  return stats;
}
