#include "col.hpp"

#include <cassert>

Col::Col() {};

void Col::reset_coloring() {
  coloring.clear();
  classes.clear();
  colorP.clear();
  colorQ.clear();
}

void Col::set_color(VertexId v, size_t pi, size_t qj, Color k) {
  assert(!is_colored(v));
  coloring[v] = k;
  classes[k].insert(v);
  colorP[pi] = k;
  if (colorQ.contains(qj)) {
    assert(colorQ[qj] == k);
  } else
    colorQ[qj] = k;
}

void Col::set_color(const DPCPInst& dpcp, VertexId vId, Color k) {
  Vertex v = vertex(vId, dpcp.get_graph());
  set_color(vId, dpcp.get_P_part(v), dpcp.get_Q_part(v), k);
}

bool Col::check_coloring(const DPCPInst& dpcp) const {
  const Graph& graph = dpcp.get_graph();

  // Return false if some P-part is uncolored
  for (auto v : boost::make_iterator_range(vertices(graph))) {
    size_t pi = dpcp.get_P_part(v);
    if (!colorP.contains(pi)) {
      std::cout << "Coloring error #2: P-part " << pi << " is uncolored"
                << std::endl;
      return false;
    }
  }

  // Return false if the coloring is not proper
  for (auto e : boost::make_iterator_range(edges(graph))) {
    Vertex u = source(e, graph);
    size_t idU = dpcp.get_current_id(u);
    Vertex v = target(e, graph);
    size_t idV = dpcp.get_current_id(v);
    if (coloring.contains(idU) && coloring.contains(idV) &&
        coloring.at(idU) == coloring.at(idV)) {
      std::cout << "Coloring error #3: " << u << " (" << dpcp.get_P_part(u)
                << "," << dpcp.get_Q_part(u) << ") and " << v << " ("
                << dpcp.get_P_part(v) << "," << dpcp.get_Q_part(v)
                << ") are adjecent and both have color " << coloring.at(idU)
                << std::endl;
      return false;
    }
  }

  return true;
}

StableEnv Col::get_stable(const DPCPInst& dpcp, Color k) const {
  StableEnv stab;
  stab.stable.reserve(classes.at(k).size());
  for (VertexId idV : classes.at(k)) {
    Vertex v = vertex(idV, dpcp.get_graph());
    stab.stable.push_back(v);
    stab.ps.insert(dpcp.get_P_part(v));
    stab.qs.insert(dpcp.get_Q_part(v));
  }
  return stab;
}

Col Col::translate_coloring(const DPCPInst& currentDpcp,
                            const DPCPInst& originalDPCP) const {
  Col dstCol;
  for (auto& [idv, k] : coloring) {
    Vertex v = vertex(idv, currentDpcp.get_graph());
    size_t originalId = currentDpcp.get_original_id(v);
    dstCol.set_color(originalDPCP, originalId, k);
  }
  return dstCol;
}

void Col::color_isolated_vertices(
    const DPCPInst& dpcp, const std::list<IsolatedVertex>& isolatedVertices) {
  for (const auto& iso : isolatedVertices) {
    VertexId id = iso.id;
    Vertex v = vertex(id, dpcp.get_graph());
    size_t pi = dpcp.get_P_part(v);
    size_t qj = dpcp.get_Q_part(v);
    // Decide color
    Color k;
    if (is_colored_Q(qj))
      k = get_color_Q(qj);
    else
      k = 0;  // First color;
    set_color(id, pi, qj, k);
  }
}

void Col::write_coloring(std::ostream& out) const {
  out << coloring.size() << " " << classes.size() << "\n";
  for (auto [idv, k] : coloring) {
    out << idv << " " << k << "\n";
  }
}
