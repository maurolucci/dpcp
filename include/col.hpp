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

  [[nodiscard]] inline const ColorClass &get_color_classes() const {
    return classes;
  };

  [[nodiscard]] inline size_t get_n_colors() const { return classes.size(); };

  [[nodiscard]] inline bool is_colored_B(TypeB b) const {
    return colorB.contains(b);
  };

  [[nodiscard]] inline bool get_color_B(TypeB b) const { return colorB[b]; };

  void reset_coloring();

  void set_color(const Vertex v, const Color k);

  [[nodiscard]] bool check_coloring(const Graph &graph) const;

  [[nodiscard]] StableEnv get_stable(const Graph &graph, const Color k) const;

private:
  Graph &graph;
  Coloring coloring;
  ColorClass classes;
  std::map<TypeA, Color> colorA;
  std::map<TypeB, Color> colorB;
};

#endif // _COL_HPP_