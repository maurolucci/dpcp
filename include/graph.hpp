#ifndef _GRAPH_HPP_
#define _GRAPH_HPP_

#include <iostream>

#include "params.hpp"

#include "boost/graph/adjacency_list.hpp"

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
using ASet = std::set<TypeA>;

using GCPGraph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, TypeB>;

std::tuple<Graph, size_t, size_t> read_dpcp_instance(std::istream &graph,
                                                     std::istream &partA,
                                                     std::istream &partB);

void get_gcp_graph(Graph &src, GCPGraph &dst, std::map<TypeB, size_t> &tyB2idB,
                   std::vector<TypeB> &idB2TyB);

void graph_copy(const Graph &src, const std::map<Vertex, size_t> &getId,
                Graph &dst);
void graph_copy(const Graph &src, Graph &dst);

void vertex_branching1(Graph &graph, Vertex v);
void vertex_branching2(Graph &graph, Vertex v);

class GraphEnv {

public:
  Graph *graphPtr;                // Graph G = (V,E) with V c AxB
  Graph &graph;                   // Graph G = (V,E) with V c AxB
  std::map<Vertex, size_t> getId; // Map from V to {0,..,|V|-1}
  size_t nA, nB;                  // |A| and |B|

  std::map<TypeA, std::vector<Vertex>> Va; // V_a = {v : v = (a,b) for some b}
  std::map<TypeB, std::vector<Vertex>> Vb; // V^b = {v : v = (a,b) for some a}

  std::map<TypeA, size_t> tyA2idA;      // Map from TypeA to idA
  std::map<TypeB, size_t> tyB2idB;      // Map from TypeB to idB
  std::vector<TypeA> idA2TyA;           // Map from idA to TypeA
  std::vector<TypeB> idB2TyB;           // Map from idB to TypeB
  std::vector<std::vector<Vertex>> snd; // V_a = {v : v = (a,b) for some b}

  bool isRoot;                    // Is the root node?
  bool isGCP;                     // Is a GCP instance?
  bool isInfeasible;              // Is the instance infeasible?
  bool hasTrivialSolution;        // Does the instance have a trivial solution?
  std::list<VertexInfo> isolated; // List of isolated vertices
  double density;                 // Density of the graph (before preprocessing)

  GraphEnv(Graph *graph, bool preprocess1 = true, bool preprocess2 = true,
           bool preprocess3 = true, bool preprocess4 = true,
           bool isRoot = false);
  ~GraphEnv();

private:
  // Intialization
  void init_graphenv();

  // Preprocessing
  void preprocess(bool preprocess1, bool preprocess2, bool preprocess3,
                  bool preprocess4);
  void init_preprocess();
  void preprocess_step1();
  void preprocess_step2();
  void preprocess_step3();
  void preprocess_step4();
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
  void add_vertex(const Vertex v, const TypeA a, const TypeB b);
};

using Pool = std::vector<StableEnv>;

#endif //_GRAPH_HPP_
