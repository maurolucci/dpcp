#ifndef _GRAPH_HPP_
#define _GRAPH_HPP_

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

void get_conflict_graph(const HGraph &hg, Graph &graph);

#endif //_GRAPH_HPP_