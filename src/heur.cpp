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

using ClockType = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::_V2::system_clock::time_point;

// Struct with the vertex information needed by the 1-step heuristic for DPCP
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
// (ii) v has a greater uncolored degree
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

// Function to select the next vertex to color using the vertex criterion VC1
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
// Criterion CC2: between all free colors, assign the one that invalidates the
// lowest number of vertices
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

    // 4. Remove every neighbor (a',b') of u with a != a' and b != b' such that
    // some vertex of Vb' has color i
    // Also, add i as an adjacent color
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
