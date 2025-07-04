#ifndef _GRAPH_HPP_
#define _GRAPH_HPP_

#include <iostream>

#include "col.hpp"
#include "params.hpp"

#include "boost/graph/adjacency_list.hpp"
#include "hglib.h"

using TypeA = size_t;
using TypeB = size_t;
struct VertexInfo {
  TypeA first;
  TypeB second;
  size_t id;
};
using Graph = boost::adjacency_list<boost::listS, boost::listS,
                                    boost::undirectedS, VertexInfo>;
using Vertex = Graph::vertex_descriptor;

using VertexVector = std::vector<Vertex>;
using VertexSet = std::set<Vertex>;

using HGraph = hglib::UndirectedHypergraph<>;
using HVertex = hglib::VertexIdType;
using HEdge = hglib::HyperedgeIdType;

using GCPGraph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, TypeB>;

void read_hypergrah(HGraph &hg, std::istream &input);

void get_conflict_graph(const HGraph &src, Graph &dst);

void get_gcp_graph(Graph &src, GCPGraph &dst, std::map<TypeB, size_t> &tyB2idB,
                   std::vector<TypeB> &idB2TyB);

void vertex_branching1(Graph &graph, Vertex v);
void vertex_branching2(Graph &graph, Vertex v);

class GraphEnv {

public:
  Graph graph;                          // Graph G = (V,E) with V c AxB
  Params &params;                       // Parameters
  size_t nA, nB;                        // |A| and |B|
  std::map<TypeA, size_t> tyA2idA;      // Map from TypeA to idA
  std::map<TypeB, size_t> tyB2idB;      // Map from TypeB to idB
  std::vector<TypeA> idA2TyA;           // Map from idA to TypeA
  std::vector<TypeB> idB2TyB;           // Map from idB to TypeB
  std::vector<std::vector<Vertex>> snd; // V_a = {v : v = (a,b) for some b}
  std::vector<std::vector<Vertex>> fst; // V^b = {v : v = (a,b) for some a}
  bool isGCP;                           // Is a GCP instance?
  bool isInfeasible;                    // Is the instance infeasible?
  std::list<VertexInfo> isolated;       // List of isolated vertices

  GraphEnv(const Graph &graph, Params &param);
  GraphEnv(const Graph &&graph, Params &param);
  ~GraphEnv();

  void color_isolated(Col &col);

private:
  // Intialization
  void init_graphenv();

  // Preprocessing
  void preprocess();
  void init_preprocess();
  void preprocess_step1();
  void preprocess_step2();
  void preprocess_step3();
};

struct StableEnv {
  VertexVector stable; // stable set
  std::set<TypeA> as;  // set of A-values from the stable set, i.e. {a \in
                       // A: (a,b) \in S for some b \in B}
  std::set<TypeB> bs;  // set of B-values from the stable set, i.e. {b \in
                       // B: (a,b) \in S for some a \in A}
  double cost;
  StableEnv();
  StableEnv(VertexVector &stable, std::set<TypeA> &as, std::set<TypeB> &bs,
            double cost);
  StableEnv(VertexVector &&stable, std::set<TypeA> &&as, std::set<TypeB> &&bs,
            double cost);

  bool check(const Graph &graph);
};

using Pool = std::vector<StableEnv>;

#endif //_GRAPH_HPP_