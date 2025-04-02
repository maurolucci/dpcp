#include "col.hpp"

Col::Col(Graph &graph) : graph(graph){};

void Col::reset_coloring() {
  coloring.clear();
  classes.clear();
}

void Col::set_color(const TypeA a, const TypeB b, const Color k) {
  // Find vertex u in the original graph such that graph[u] = (a,b)
  Vertex u = num_vertices(graph);
  for (auto v : boost::make_iterator_range(vertices(graph)))
    if (graph[v].first == a && graph[v].second == b) {
      u = v;
      break;
    }
  assert(u < num_vertices(graph));
  coloring[u] = k;
  classes[k].insert(u);
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
        coloring.at(u) == coloring.at(v))
      return false;
  }

  return true;
}