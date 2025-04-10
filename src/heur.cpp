#include "heur.hpp"
#include "random.hpp"
extern "C" {
#include "color.h"
}

#include <chrono>
#include <limits>
#include <random>

using ClockType = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::_V2::system_clock::time_point;

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
  std::vector<VertexVector> adj;        // Adjacent lists of the subgraph
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
      adj.push_back(VertexVector());
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
        col.set_color(v, k);
    }
  // assert(col.check_coloring());
  return;
}

Stats heur_solve(const GraphEnv &genv, const std::vector<TypeA> &as, Col &col,
                 size_t repetitions) {
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
  }
  assert(col.check_coloring(genv.graph));
  TimePoint end = ClockType::now();
  return Stats{0,
               0,
               FEASIBLE,
               std::chrono::duration<double>(end - start).count(),
               0,
               std::numeric_limits<double>::lowest(),
               static_cast<double>(col.get_n_colors()),
               std::numeric_limits<double>::max()};
}
