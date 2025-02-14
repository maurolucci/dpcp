#ifndef _GRAPH_HPP_
#define _GRAPH_HPP_

#include "boost/graph/adjacency_list.hpp"
#include "hglib.h"

using HyperGraph = hglib::UndirectedHypergraph<>;
using Vertex = hglib::VertexIdType;
using Hyperedge = hglib::HyperedgeIdType;
using VertexSet = std::set<Vertex>;

using PSet = std::pair<Hyperedge, Vertex>;
using ConflictGraph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, PSet>;
using CGVertex = ConflictGraph::vertex_descriptor;

void init_conflict_graph(const HyperGraph &hg, ConflictGraph &cg);

#endif //_GRAPH_HPP_