#include "col.hpp"

Col::Col(){};

void Col::reset_coloring() {
  coloring.clear();
  classes.clear();
  colorA.clear();
  colorB.clear();
}

void Col::set_color(const Graph &graph, const Vertex v, const Color k) {
  TypeA a = graph[v].first;
  TypeB b = graph[v].second;
  coloring[v] = k;
  classes[k].insert(v);
  colorA[a] = k;
  if (colorB.contains(b))
    // Al the vertices of V^b must have the same color
    assert(colorB[b] == k);
  else {
    colorB[b] = k;
  }
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
      std::cout << "Coloring error #3: " << u << " (" << graph[u].first << ","
                << graph[u].second << ") and " << v << " (" << graph[v].first
                << "," << graph[v].second
                << ") are adjecent and both have color " << coloring.at(u)
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

void Col::translate_coloring(const Graph &srcGraph, const Graph &dstGraph,
                             Col &dstCol) {
  for (auto &[v, k] : coloring) {
    size_t id = srcGraph[v].id;
    // Find original vertex
    auto u = vertex(id, dstGraph);
    dstCol.set_color(dstGraph, u, k);
  }
}

void Col::color_isolated_vertices(std::list<VertexInfo> &isolated, Col &dstCol,
                                  const Graph &dstGraph) {
  for (auto [a, b, id] : isolated) {
    // Find original vertex
    auto v = vertex(id, dstGraph);
    // Decide color
    Color k;
    if (dstCol.is_colored_B(b))
      k = dstCol.get_color_B(b);
    else
      k = 0; // First color;
    dstCol.set_color(dstGraph, v, k);
  }
}

void Col::write_coloring(std::ostream &out) const {
  out << coloring.size() << " " << classes.size() << "\n";
  for (auto [v, k] : coloring) {
    out << v << " " << k << "\n";
  }
}