#include "feas.hpp"
#include "heur.hpp"

extern "C" {
#include "color.h"
#include "color_private.h"
#include "mwis.h"
}
#include <cfloat>
#include <chrono>
#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>

#define FEASIBILITY_EPSILON 0.00001 // 10e-5

Stats dpcp_decide_feasibility_enumerative(const GraphEnv &_genv, Col &col,
                                          std::ostream &log) {

  // Initial time instant
  auto start = std::chrono::high_resolution_clock::now();

  // First make a copy of the graph
  Graph g;
  graph_copy(_genv.graph, g);
  GraphEnv genv(&g, true, false, false, false, true);

  // Remove edges whose endpoints do not belongs to the same Va and Vb
  Graph::edge_iterator ei, ei_end, next;
  boost::tie(ei, ei_end) = edges(g);
  for (next = ei; ei != ei_end; ei = next) {
    ++next;
    Vertex u = source(*ei, g);
    Vertex v = target(*ei, g);
    if (genv.graph[u].first != genv.graph[v].first &&
        genv.graph[u].second != genv.graph[v].second) {
      remove_edge(*ei, g);
    }
  }

  // Initalize stable environment
  MWISenv *mwis_env = NULL;
  COLORstable_initenv(&mwis_env, NULL, 0);

  // Intialize vectors of weights
  COLORNWT *mwis_pi = NULL;
  mwis_pi = (COLORNWT *)COLOR_SAFE_MALLOC(num_vertices(g), COLORNWT);
  for (size_t i = 0; i < num_vertices(g); ++i)
    mwis_pi[i] = 1;
  COLORNWT mwis_pi_scalef = INT_MAX; // Force optimality

  // Initialize edge array
  int ecount = 0;
  int *elist = (int *)malloc(sizeof(int) * 2 * num_edges(g));
  for (auto e : boost::make_iterator_range(edges(g))) {
    elist[2 * ecount] = genv.getId[source(e, g)];
    elist[2 * ecount++ + 1] = genv.getId[target(e, g)];
  }

  // Solve the MWIS problem up to optimality
  COLORset *newsets = NULL;
  int nnewsets = 0;
  COLORstable_wrapper(&mwis_env, &newsets, &nnewsets, num_vertices(genv.graph),
                      ecount, elist, mwis_pi, mwis_pi_scalef, 0, 0, 2);

  assert(nnewsets > 0);

  // If the maximum stable set has size |A|, then we have found a feasible
  // solution of DPCP
  if (newsets[0].count == static_cast<int>(genv.nA)) {
    // First, find selected vertices
    VertexVector selected;
    std::map<TypeB, std::set<TypeB>> adj;
    for (int i = 0; i < newsets[0].count; ++i) {
      int vi = newsets[0].members[i];
      Vertex v = vertex(vi, _genv.graph);
      // Add adjacencies
      for (auto u : selected) {
        // Warning: adjacencies must be checked in the original graph, as the
        // new graph may not contain all edges
        if (edge(u, v, _genv.graph).second) {
          TypeB bb = _genv.graph[u].second;
          TypeB b = _genv.graph[v].second;
          if (b < bb)
            adj[b].insert(bb);
          else if (bb > b)
            adj[bb].insert(b);
        }
      }
      selected.push_back(v);
    }
    // Then, color them
    dpcp_dsatur_heur(_genv, selected, adj, col);
    assert(col.check_coloring(_genv.graph));
  }

  // Save stats
  Stats stats;
  stats.state =
      newsets[0].count == static_cast<int>(genv.nA) ? FEASIBLE : INFEASIBLE;
  stats.time = std::chrono::duration<double>(
                   std::chrono::high_resolution_clock::now() - start)
                   .count();
  stats.lb = newsets[0].count;

  // Free memory
  for (int i = 0; i < nnewsets; ++i)
    free(newsets[i].members);
  free(newsets);
  free(elist);
  free(mwis_pi);
  COLORstable_freeenv(&mwis_env);

  return stats;
}

// This is the class implementing the generic callback interface.
class EarlyStopCallback : public IloCplex::Callback::Function {

private:
  size_t nA;
  std::ostream &log;

public:
  // Constructor with data.
  EarlyStopCallback(size_t nA, std::ostream &log) : nA(nA), log(log){};

  // Check if we can stop early: abort if the best bound is below nA
  inline void check_early_stop(const IloCplex::Callback::Context &context) {

    if (context.inGlobalProgress())
      if (context.getDoubleInfo(IloCplex::Callback::Context::Info::BestBound) <
          nA - FEASIBILITY_EPSILON) {
        log << "Early stop: infeasibility detected" << std::endl;
        context.abort();
      }
  }

  // This is the function that we have to implement and that CPLEX will call
  // during the solution process at the places that we asked for.
  virtual void invoke(const IloCplex::Callback::Context &context) ILO_OVERRIDE {
    if (context.inGlobalProgress())
      check_early_stop(context);
  }

  /// Destructor
  ~EarlyStopCallback(){};
};

Stats dpcp_decide_feasibility_ilp(const GraphEnv &_genv, Col &col,
                                  int timeLimit, std::ostream &log) {

  // Initial time instant
  auto start = std::chrono::high_resolution_clock::now();

  // First make a copy of the graph
  Graph g;
  graph_copy(_genv.graph, g);
  GraphEnv genv(&g, true, false, false, false, true);

  // Remove edges whose endpoints do not belongs to the same Va and Vb
  Graph::edge_iterator ei, ei_end, next;
  boost::tie(ei, ei_end) = edges(g);
  for (next = ei; ei != ei_end; ei = next) {
    ++next;
    Vertex u = source(*ei, g);
    Vertex v = target(*ei, g);
    if (genv.graph[u].first != genv.graph[v].first &&
        genv.graph[u].second != genv.graph[v].second) {
      remove_edge(*ei, g);
    }
  }

  // CPLEX environment
  IloEnv cxenv;
  IloModel cxmodel(cxenv);
  IloObjective cxobj(cxenv);
  IloNumVarArray x(cxenv, num_vertices(g));
  IloConstraintArray cxcons(cxenv);
  IloCplex cplex(cxenv);

  // Variables
  for (auto v : boost::make_iterator_range(vertices(g))) {
    char name[100];
    snprintf(name, sizeof(name), "x_%ld_%ld", g[v].first, g[v].second);
    x[genv.getId[v]] = IloBoolVar(cxenv, name);
  }

  // Objective
  IloExpr obj(cxenv);
  for (auto v : boost::make_iterator_range(vertices(g)))
    obj += x[genv.getId[v]] * 1.0;
  cxobj = IloMaximize(cxenv, obj);
  cxmodel.add(cxobj);
  obj.end();

  // Constraints
  for (auto e : boost::make_iterator_range(edges(g))) {
    IloExpr expr(cxenv);
    expr += x[genv.getId[source(e, g)]] + x[genv.getId[target(e, g)]];
    cxcons.add(expr <= 1);
    expr.end();
  }
  cxmodel.add(cxcons);

  // Set parameters and extract the model
  cplex.extract(cxmodel);
  cplex.setOut(log);
  cplex.setDefaults();
  cplex.setParam(IloCplex::Param::Parallel, 1); // Deterministic mode
  cplex.setParam(IloCplex::Param::Threads, 1);  // Single thread
  cplex.setParam(IloCplex::Param::TimeLimit, timeLimit);

  // Create the callback object
  EarlyStopCallback cb(genv.nA, log);
  // Now we get to setting up the generic callback.
  // We instantiate a ThresholdCallback and set the contextMask parameter.
  CPXLONG contextMask = IloCplex::Callback::Context::Id::GlobalProgress;
  // If contextMask is not zero we add the callback.
  if (contextMask != 0)
    cplex.use(&cb, contextMask);

  // Finally, we can solve the model
  cplex.solve();

  // If optimal and obj value >= nA, then we have a feasible solution of DPCP
  // We must find the coloring in the original graph
  if (cplex.getCplexStatus() == IloCplex::Optimal &&
      static_cast<size_t>(cplex.getObjValue()) >= genv.nA) {
    // Get variable values
    IloNumArray vals(cxenv);
    cplex.getValues(vals, x);
    // First, find selected vertices
    VertexVector selected;
    std::map<TypeB, std::set<TypeB>> adj;
    for (auto vv : boost::make_iterator_range(vertices(g))) {
      Vertex v = vertex(genv.getId[vv], _genv.graph);
      if (vals[genv.getId[v]] > 0.5) {
        // Add adjacencies
        for (auto u : selected) {
          if (edge(u, v, _genv.graph).second) {
            TypeB bb = _genv.graph[u].second;
            TypeB b = _genv.graph[v].second;
            if (b < bb)
              adj[b].insert(bb);
            else if (bb > b)
              adj[bb].insert(b);
          }
        }
        // Add new selected vertex
        selected.push_back(v);
      }
    }
    // Then, color them
    dpcp_dsatur_heur(_genv, selected, adj, col);
    assert(col.check_coloring(_genv.graph));
    vals.end();
  }

  // Save stats
  Stats stats;
  if (cplex.getCplexStatus() == IloCplex::Optimal) {
    stats.state = static_cast<size_t>(cplex.getObjValue()) >= genv.nA
                      ? FEASIBLE
                      : INFEASIBLE;
    stats.lb = cplex.getObjValue();
  } else if (cplex.getCplexStatus() == IloCplex::AbortUser) {
    stats.state = INFEASIBLE;
    stats.lb = cplex.getBestObjValue();
  } else {
    stats.state = UNKNOWN;
  }
  stats.time = std::chrono::duration<double>(
                   std::chrono::high_resolution_clock::now() - start)
                   .count();
  stats.nodes = cplex.getNnodes();

  // Free memory
  cxobj.end();
  for (size_t v = 0; v < num_vertices(g); ++v)
    x[v].end();
  cxcons.end();
  cxmodel.end();
  cplex.end();
  cxenv.end();

  return stats;
}