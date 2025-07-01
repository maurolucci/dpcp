#include "col.hpp"

Col::Col(){};

void Col::reset_coloring() {
  coloring.clear();
  classes.clear();
}

// void Col::set_color(const TypeA a, const TypeB b, const Color k) {
//   // Find vertex u in the original graph such that graph[u] = (a,b)
//   Vertex u = num_vertices(graph);
//   for (auto v : boost::make_iterator_range(vertices(graph)))
//     if (graph[v].first == a && graph[v].second == b) {
//       u = v;
//       break;
//     }
//   assert(u < num_vertices(graph));
//   coloring[u] = k;
//   classes[k].insert(u);
// }

void Col::set_color(const Vertex v, const Color k) {
  coloring[v] = k;
  classes[k].insert(v);
}

bool Col::check_coloring(const Graph &graph) const {

  std::map<TypeA, Color> colorA;
  std::map<TypeB, Color> colorB;

  // Initialize the color of every element of A and B to -1
  for (auto v : boost::make_iterator_range(vertices(graph))) {
    auto [a, b] = graph[v];
    colorA[a] = -1;
    colorB[b] = -1;
  }

  // Get final colors
  // Return false if some b \in B has more than one color
  for (auto [v, k] : coloring) {
    auto [a, b] = graph[v];
    colorA[a] = k;
    if (colorB[b] != -1 && colorB[b] != k) {
      std::cout << "Coloring error #1: " << b << " has colors " << colorB[b]
                << " and " << k << std::endl;
      return false;
    }
    colorB[b] = k;
  }

  // Return false if some a \in A is uncolored
  for (auto [a, k] : colorA)
    if (k == -1) {
      std::cout << "Coloring error #2: " << a << " is uncolored" << std::endl;
      return false;
    }

  // Return false if the coloring is not proper
  for (auto e : boost::make_iterator_range(edges(graph))) {
    auto u = source(e, graph);
    auto v = target(e, graph);
    if (coloring.contains(u) && coloring.contains(v) &&
        coloring.at(u) == coloring.at(v)) {
      std::cout << "Coloring error #3: " << u << " and " << v
                << " are adjecent and both have color " << coloring.at(u)
                << std::endl;
      return false;
    }
  }

  return true;
}

StableEnv Col::get_stable(const Graph &graph, const Color k) const {
  StableEnv stab;
  stab.stable.reserve(classes.at(k).size());
  for (Vertex v : classes.at(k)) {
    stab.stable.push_back(v);
    stab.as.insert(graph[v].first);
    stab.bs.insert(graph[v].second);
  }
  return stab;
}