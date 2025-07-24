#include "heur.hpp"
#include "graph.hpp"
#include "random.hpp"
extern "C" {
#include "color.h"
}

#include <chrono>
#include <limits>
#include <random>

using ClockType = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::_V2::system_clock::time_point;

class CustomComparator {
public:
  CustomComparator(const GraphEnv &genv, const std::vector<bool> &state,
                   const std::vector<int> &nV)
      : genv(genv), state(state), nV(nV) {}

  size_t degree(Vertex v) const {
    size_t ret = 0;
    for (Vertex w :
         boost::make_iterator_range(adjacent_vertices(v, genv.graph)))
      if (state[genv.getId.at(v)])
        ret++;
    return ret;
  }

  Vertex select() {
    auto x = std::min_element(vertices(genv.graph).first,
                              vertices(genv.graph).second);
    return *x;
  }

  bool operator()(Vertex v, Vertex u) const {

    // Get the a component of each vertex
    TypeA av = genv.graph[v].first, au = genv.graph[u].first;

    // First criterion: smallest |Va|
    if (nV[genv.tyA2idA.at(av)] < nV[genv.tyA2idA.at(au)])
      return true;
    else if (nV[genv.tyA2idA.at(av)] > nV[genv.tyA2idA.at(au)])
      return false;

    // Second criterion: smallest |N(w)|
    size_t degreev = degree(v);
    size_t degreeu = degree(u);
    return degreev <= degreeu;
  }

private:
  const GraphEnv &genv;
  const std::vector<bool> &state;
  const std::vector<int> &nV;
};

Vertex select_vertex(const GraphEnv &genv, const std::vector<bool> &state) {}

// One-step heuristic for DPCP
void dpcp_heur_1_step(const GraphEnv &genv, Col &col, VertexSelector getVertex,
                      StableSetConstructor buildStab) {

  // Map from vertex to state. Possible states are:
  //  * false: unavailable (colored or removed)
  //  * true: available (uncolored)
  std::vector<bool> state(num_vertices(genv.graph), true);

  // First color
  Color k = 0;

  // Map from a to the cardinality of V_a
  std::vector<size_t> nV(genv.nA);
  for (size_t iA = 0; iA < genv.nA; ++iA)
    nV[iA] = genv.snd[iA].size();

  // Number of vertices to finish
  size_t n_rest = genv.nA;

  while (n_rest > 0) {

    // Chose a vertex
    Vertex v = getVertex(genv, state);

    // Build an stable set containing v
    VertexVector stable = buildStab(genv, state, v);

    // Iterate over the vertices of the stable set
    for (Vertex u : stable) {

      // Color u with k
      col.set_color(genv.graph, u, k);

      // Get the a and b components of u
      TypeA a = genv.graph[u].first;
      TypeB b = genv.graph[u].second;

      // Remove vertices of Va
      for (Vertex w : genv.snd.at(genv.tyA2idA.at(a)))
        state[genv.getId.at(w)] = false;
      nV[genv.tyA2idA.at(a)] = std::numeric_limits<size_t>::max();
      n_rest--;

      // Iterate over the neighbors of u
      for (Vertex w :
           boost::make_iterator_range(adjacent_vertices(u, genv.graph))) {

        // Remove w if it also has b
        if (!state[genv.getId.at(w)] || genv.graph[w].second != b)
          continue;
        state[genv.getId.at(w)] = false;

        // Get the a component of w
        TypeA aw = genv.graph[w].first;

        // If aw has no more candidates, report failure
        nV[genv.tyA2idA.at(aw)]--;
        if (nV[genv.tyA2idA.at(aw)] > 0)
          continue;
        col.reset_coloring();
        return;
      }
    }

    // Next color
    ++k;
  }

  assert(col.check_coloring(genv.graph));
  return;
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
