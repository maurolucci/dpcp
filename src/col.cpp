#include "col.hpp"

void Col::set_color(const Vertex v, const Color k) {
  if (coloring.contains(v))
    classes[coloring[v]].erase(v);
  coloring[v] = k;
  classes[k].insert(v);
}

bool Col::check_coloring() const {

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
    if (colorB[b] != -1 && colorB[b] != k)
      return false;
    colorB[b] = k;
  }

  // Return false if some a \in A is uncolored
  for (auto [a, k] : colorA)
    if (k == -1)
      return false;

  // Return false if the coloring is not proper
  for (auto e : boost::make_iterator_range(edges(graph))) {
    auto u = source(e, graph);
    auto v = target(e, graph);
    if (coloring.contains(u) && coloring.contains(v) &&
        coloring[u] != coloring[v])
      return false;
  }

  return true;
}