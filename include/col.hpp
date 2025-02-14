#ifndef _COL_HPP_
#define _COL_HPP_

#include "graph.hpp"

#include <vector>

using Color = int;
using Coloring = std::map<Vertex, Color>;
using VertexSet = std::set<Vertex>;
using ColorClass = std::map<Color, VertexSet>;

class Col {

public:
  [[nodiscard]] inline const Coloring &get_coloring() const {
    return coloring;
  };

  [[nodiscard]] inline size_t get_n_colors() const { return classes.size(); };

  void set_color(const Vertex v, const Color k);

  [[nodiscard]] bool check_coloring() const;

private:
  Graph &graph;
  Coloring coloring;
  ColorClass classes;
};

#endif // _COL_HPP_