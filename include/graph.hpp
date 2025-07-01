#ifndef _GRAPH_HPP_
#define _GRAPH_HPP_

#include <iostream>

#include "boost/graph/adjacency_list.hpp"
#include "hglib.h"

using TypeA = size_t;
using TypeB = size_t;
using Pair = std::pair<TypeA, TypeB>;
using Graph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Pair>;
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

struct GraphEnv {
  Graph graph;                     // Graph G = (V,E) with V c AxB
  size_t nA, nB;                   // |A| and |B|
  std::map<TypeA, size_t> tyA2idA; // Map from TypeA to idA
  std::map<TypeB, size_t> tyB2idB; // Map from TypeB to idB
  std::vector<TypeA> idA2TyA;      // Map from idA to TypeA
  std::vector<TypeB> idB2TyB;      // Map from idB to TypeB
  std::vector<std::vector<Vertex>>
      snd; // Map from idA to subset of Vertex:
           // snd[i_a] = {(a,b) in V: tyA2idA[a] = i_a}
  std::vector<std::vector<Vertex>>
      fst;    // Map from idB to subset of Vertex:
              // fst[i_b] = {(a,b) in V: tyB2idB[b] = i_b}
  bool isGCP; // Whether the instance is a graph coloring instance, i.e.
  // |snd[a]| = 1 forall a \in A

  GraphEnv(const Graph &graph);
  GraphEnv(const Graph &&graph);
  ~GraphEnv();
};

struct StableEnv {
  VertexVector stable; // stable set
  std::set<TypeA> as;  // set of A-values from the stable set, i.e. {a \in
                       // A: (a,b) \in S for some b \in B}
  std::set<TypeB> bs;  // set of B-values from the stable set, i.e. {b \in
                       // B: (a,b) \in S for some a \in A}
  double cost;
  StableEnv();
  StableEnv(VertexVector &stable, std::set<TypeB> &as, std::set<TypeB> &bs,
            double cost);
  StableEnv(VertexVector &&stable, std::set<TypeB> &&as, std::set<TypeB> &&bs,
            double cost);

  bool check(const Graph &graph);
};

using Pool = std::vector<StableEnv>;

#endif //_GRAPH_HPP_