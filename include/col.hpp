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
  Col();

  [[nodiscard]] inline const Coloring &get_coloring() const {
    return coloring;
  };

  [[nodiscard]] inline const ColorClass &get_color_classes() const {
    return classes;
  };

  [[nodiscard]] inline size_t get_n_colors() const { return classes.size(); };

  [[nodiscard]] inline Color get_color(const Vertex v) const {
    return coloring.at(v);
  };

  void reset_coloring();

  void set_color(const Vertex v, const Color k);

  [[nodiscard]] bool check_coloring(const Graph &graph) const;

private:
  Coloring coloring;
  ColorClass classes;
};

#endif // _COL_HPP_