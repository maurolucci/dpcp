#ifndef _COL_HPP_
#define _COL_HPP_

#include "graph.hpp"

#include <map>
#include <vector>

using Color = int;
using Coloring = std::map<Vertex, Color>;
using VertexSet = std::set<Vertex>;
using ColorClass = std::map<Color, VertexSet>;

class Col {

public:
  Col(Graph &graph);

  [[nodiscard]] inline const Coloring &get_coloring() const {
    return coloring;
  };

  [[nodiscard]] inline size_t get_n_colors() const { return classes.size(); };

  void reset_coloring();

  void set_color(const TypeA a, const TypeB b, const Color k);

  [[nodiscard]] bool check_coloring() const;

private:
  const Graph graph;
  Coloring coloring;
  ColorClass classes;
};

#endif // _COL_HPP_