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

  [[nodiscard]] inline bool is_colored_B(TypeB b) const {
    return colorB.contains(b);
  };

  [[nodiscard]] inline Color get_color_B(TypeB b) const {
    return colorB.at(b);
  };

  void reset_coloring();

  void set_color(const Graph &graph, const Vertex v, const Color k);

  [[nodiscard]] bool check_coloring(const Graph &graph) const;

  [[nodiscard]] StableEnv get_stable(const Graph &graph, const Color k) const;

  void translate_coloring(const Graph &srcGraph, const Graph &dstGraph,
                          Col &dstCol);

  void color_isolated_vertices(std::list<VertexInfo> &isolated, Col &dstCol,
                               const Graph &dstGraph);

  void write_coloring(std::ostream &out) const;

private:
  Coloring coloring;
  ColorClass classes;
  std::map<TypeA, Color> colorA;
  std::map<TypeB, Color> colorB;
};

#endif // _COL_HPP_