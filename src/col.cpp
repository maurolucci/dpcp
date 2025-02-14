#include "col.hpp"

void Col::set_color(const Vertex v, const Color k) {
  if (coloring.contains(v))
    classes[coloring[v]].erase(v);
  coloring[v] = k;
  classes[k].insert(v);
}

bool Col::check_coloring() const {
  for (const auto &e : hg.hyperedges()) {
    std::map<Color, size_t> counter;
    for (const auto &v : *hg.impliedVertices(e->id()))
      if (coloring.contains(v->id()))
        counter[coloring.at(v->id())]++;
    auto it = std::find_if(
        counter.begin(), counter.end(),
        [](const std::pair<Color, size_t> &p) { return p.second == 1; });
    if (it == counter.end())
      return false;
  }
  return true;
}