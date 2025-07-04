#include "col.hpp"

Col::Col(Graph &graph) : graph(graph) {};

void Col::reset_coloring() {
  coloring.clear();
  classes.clear();
  colorA.clear();
  colorB.clear();
}

void Col::set_color(const Vertex v, const Color k) {
  TypeA a = graph[v].first;
  TypeB b = graph[v].second;
  coloring[v] = k;
  classes[k].insert(v);
  colorA[a] = k;
  if (colorB.contains(b))
    assert(colorB[b] == k);
  else
    colorB[b] = k;
}

bool Col::check_coloring(const Graph &graph) const {

  // Return false if some b \in B has more than one color
  // Already checked by set_color

  // Return false if some a \in A is uncolored
  for (auto v : boost::make_iterator_range(vertices(graph))) {
    TypeA a = graph[v].first;
    if (!colorA.contains(a)) {
      std::cout << "Coloring error #2: " << a << " is uncolored" << std::endl;
      return false;
    }
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