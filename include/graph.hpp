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

using HGraph = hglib::UndirectedHypergraph<>;
using HVertex = hglib::VertexIdType;
using HEdge = hglib::HyperedgeIdType;

using GCPGraph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, TypeB>;

void read_hypergrah(HGraph &hg, std::istream &input);

void get_conflict_graph(const HGraph &src, Graph &dst);

void get_gcp_graph(Graph &src, GCPGraph &dst, std::map<TypeB, size_t> &tyB2idB,
                   std::vector<TypeB> &idB2TyB);

#endif //_GRAPH_HPP_