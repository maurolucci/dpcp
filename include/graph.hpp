#ifndef _GRAPH_HPP_
#define _GRAPH_HPP_

#include <iostream>

#include "boost/graph/adjacency_list.hpp"
#include "params.hpp"

struct VertexInfo {
  size_t id;  // Original vertex id in {0, ..., |V|-1}
};

struct IsolatedVertex {
  size_t id;  // Original vertex id
};
using Graph = boost::adjacency_list<boost::listS, boost::listS,
                                    boost::undirectedS, VertexInfo>;
using Vertex = Graph::vertex_descriptor;

using VertexVector = std::vector<Vertex>;
using VertexSet = std::set<Vertex>;

using GCPGraph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, size_t>;

class DPCPInst;

template <typename T>
using VertexMap = std::map<Vertex, T>;
using Partition = std::vector<VertexVector>;

void graph_copy(const Graph& src, const VertexMap<size_t>& vertex2CurrentId,
                Graph& dst);
void graph_copy(const Graph& src, Graph& dst);

void read_dpcp_instance(std::istream& graphStream, std::istream& partP,
                        std::istream& partQ, Graph& graph, Partition& P,
                        Partition& Q);

class DPCPInst {
 public:
  // Constructors and destructor
  DPCPInst(const Graph& graph, const Partition& P, const Partition& Q);
  DPCPInst(const DPCPInst& dpcp);

  // Move constructor and move assignment operator are deleted to avoid
  // accidental moves that can lead to dangling references in the graph and
  // partitions.
  DPCPInst(DPCPInst&& dpcp) = delete;
  DPCPInst& operator=(DPCPInst&& dpcp) = delete;

  ~DPCPInst();

  // Check the consistency of the instance
  bool check_consistency() const;

  // Getters
  [[nodiscard]] Graph& get_graph() { return graph; }
  [[nodiscard]] const Graph& get_graph() const { return graph; }
  [[nodiscard]] const VertexMap<size_t>& get_vertex2CurrentId() const {
    return vertex2CurrentId;
  }
  [[nodiscard]] const VertexMap<size_t>& get_vertex2Ppart() const {
    return vertex2Ppart;
  }
  [[nodiscard]] const VertexMap<size_t>& get_vertex2Qpart() const {
    return vertex2Qpart;
  }
  [[nodiscard]] const Partition& get_P() const { return P; }
  [[nodiscard]] const Partition& get_Q() const { return Q; }
  [[nodiscard]] const std::list<IsolatedVertex>& get_isolated_vertices() const {
    return isolated;
  }

  // Other methods
  [[nodiscard]] size_t get_nP() const { return P.size(); }
  [[nodiscard]] size_t get_nQ() const { return Q.size(); }

  [[nodiscard]] bool is_gcp_instance() const { return isGCP; }
  [[nodiscard]] bool is_infeasible_instance() const { return isInfeasible; }
  [[nodiscard]] bool has_trivial_solution() const { return hasTrivialSolution; }
  [[nodiscard]] double get_density() const { return density; }

  // Get the P-part, Q-part, current id, and original id of a vertex
  [[nodiscard]] size_t get_P_part(Vertex v) const { return vertex2Ppart.at(v); }
  [[nodiscard]] size_t get_Q_part(Vertex v) const { return vertex2Qpart.at(v); }
  [[nodiscard]] size_t get_current_id(Vertex v) const {
    return vertex2CurrentId.at(v);
  }
  [[nodiscard]] size_t get_original_id(Vertex v) const { return graph[v].id; }
  [[nodiscard]] bool has_vertex(Vertex v) const {
    return vertex2CurrentId.contains(v);
  }

  // Get the GCP instance corresponding to the current DPCP instance
  [[nodiscard]] GCPGraph get_gcp_graph() const;

  // Remove one vertex and keep all internal maps/partitions consistent.
  void remove_vertex(Vertex v);

  // Preprocessing
  void preprocess(bool clique = false);

  // Branching decisions
  void preselect_vertex(Vertex v);
  void forbid_vertex(Vertex v);

 private:
  Graph graph;                         // Graph G = (V,E)
  VertexMap<size_t> vertex2CurrentId;  // Map from V to {0,..,|V|-1}
  Partition P;  // P: partition of V into P-parts, P[pi] = vertices in P-part pi
  Partition Q;  // Q: partition of V into Q-parts, Q[qj] = vertices in Q-part qj
  VertexMap<size_t> vertex2Ppart;      // Map from V to P-part index
  VertexMap<size_t> vertex2Qpart;      // Map from V to Q-part index
  std::list<IsolatedVertex> isolated;  // List of isolated vertices

  bool isGCP;               // Is a GCP instance?
  bool isInfeasible;        // Is the instance infeasible?
  bool hasTrivialSolution;  // Does the instance have a trivial solution?
  double density;           // Density of the graph (before preprocessing)

  // Internal preprocessing steps
  void preprocess_step1();
  void preprocess_step2();
  void preprocess_step3();
  void preprocess_step4();
};

struct StableEnv {
  VertexVector stable;  // stable set
  std::set<size_t> ps;  // set of P-parts from the stable set
  std::set<size_t> qs;  // set of Q-parts from the stable set
  double cost;
  StableEnv();
  StableEnv(VertexVector& stable, std::set<size_t>& ps, std::set<size_t>& qs,
            double cost);
  StableEnv(VertexVector&& stable, std::set<size_t>&& ps, std::set<size_t>&& qs,
            double cost);

  bool check(const Graph& graph);
  void add_vertex(const Vertex v, size_t pi, size_t qj);
};

#endif  //_GRAPH_HPP_
