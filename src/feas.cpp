#include "feas.hpp"

#include "heur.hpp"

extern "C" {
#include "color.h"
#include "color_private.h"
#include "mwis.h"
}
#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>

#include <cfloat>
#include <chrono>

#define FEASIBILITY_EPSILON 0.00001  // 10e-5

Stats dpcp_decide_feasibility_enumerative(DPCPInst& inputDpcp, Col& col,
                                          std::ostream& log) {
  // Initial time instant
  auto start = std::chrono::high_resolution_clock::now();

  // Remove edges whose endpoints do not belongs to the same Pi and Qj
  Graph::edge_iterator ei, ei_end, next;
  boost::tie(ei, ei_end) = edges(g);
  for (next = ei; ei != ei_end; ei = next) {
    ++next;
    Vertex u = source(*ei, g);
    Vertex v = target(*ei, g);
    if (dpcp.get_P_part(u) != dpcp.get_P_part(v) &&
        dpcp.get_Q_part(u) != dpcp.get_Q_part(v)) {
      remove_edge(*ei, g);
    }
  }

  // Initalize stable environment
  MWISenv* mwis_env = NULL;
  COLORstable_initenv(&mwis_env, NULL, 0);

  // Intialize vectors of weights
  COLORNWT* mwis_pi = NULL;
  mwis_pi = (COLORNWT*)COLOR_SAFE_MALLOC(num_vertices(g), COLORNWT);
  for (size_t i = 0; i < num_vertices(g); ++i) mwis_pi[i] = 1;
  COLORNWT mwis_pi_scalef = INT_MAX;  // Force optimality

  // Initialize edge array
  int ecount = 0;
  int* elist = (int*)malloc(sizeof(int) * 2 * num_edges(g));
  for (auto e : boost::make_iterator_range(edges(g))) {
    elist[2 * ecount] = dpcp.get_current_id(source(e, g));
    elist[2 * ecount++ + 1] = dpcp.get_current_id(target(e, g));
  }

  // Solve the MWIS problem up to optimality
  COLORset* newsets = NULL;
  int nnewsets = 0;
  COLORstable_wrapper(&mwis_env, &newsets, &nnewsets,
                      num_vertices(dpcp.get_graph()), ecount, elist, mwis_pi,
                      mwis_pi_scalef, 0, 0, 2);

  assert(nnewsets > 0);

  // If the maximum stable set has size |A|, then we have found a feasible
  // solution of DPCP
  if (newsets[0].count == static_cast<int>(dpcp.get_nP())) {
    // First, find selected vertices
    VertexVector selected;
    std::map<size_t, std::set<size_t>> adj;
    for (int i = 0; i < newsets[0].count; ++i) {
      int vi = newsets[0].members[i];
      Vertex v = vertex(vi, inputDpcp.get_graph());
      // Add adjacencies
      for (auto u : selected) {
        // Warning: adjacencies must be checked in the original graph, as the
        // new graph may not contain all edges
        if (edge(u, v, inputDpcp.get_graph()).second) {
          size_t qj1 = inputDpcp.get_Q_part(u);
          size_t qj2 = inputDpcp.get_Q_part(v);
          if (qj1 < qj2)
            adj[qj1].insert(qj2);
          else if (qj2 > qj1)
            adj[qj2].insert(qj1);
        }
      }
      selected.push_back(v);
    }
    // Then, color them
    dpcp_dsatur_heur(inputDpcp, selected, adj, col);
    assert(col.check_coloring(inputDpcp));
  }

  // Save stats
  Stats stats;
  stats.state = newsets[0].count == static_cast<int>(dpcp.get_nP())
                    ? FEASIBLE
                    : INFEASIBLE;
  stats.time = std::chrono::duration<double>(
                   std::chrono::high_resolution_clock::now() - start)
                   .count();
  stats.lb = newsets[0].count;

  // Free memory
  for (int i = 0; i < nnewsets; ++i) free(newsets[i].members);
  free(newsets);
  free(elist);
  free(mwis_pi);
  COLORstable_freeenv(&mwis_env);

  return stats;
}

// This is the class implementing the generic callback interface.
class EarlyStopCallback : public IloCplex::Callback::Function {
 private:
  size_t nP;
  std::ostream& log;

 public:
  // Constructor with data.
  EarlyStopCallback(size_t nP, std::ostream& log) : nP(nP), log(log) {};

  // Check if we can stop early: abort if the best bound is below nP
  inline void check_early_stop(const IloCplex::Callback::Context& context) {
    if (context.inGlobalProgress())
      if (context.getDoubleInfo(IloCplex::Callback::Context::Info::BestBound) <
          nP - FEASIBILITY_EPSILON) {
        log << "Early stop: infeasibility detected" << std::endl;
        context.abort();
      }
  }

  // This is the function that we have to implement and that CPLEX will call
  // during the solution process at the places that we asked for.
  virtual void invoke(const IloCplex::Callback::Context& context) ILO_OVERRIDE {
    if (context.inGlobalProgress()) check_early_stop(context);
  }

  /// Destructor
  ~EarlyStopCallback() {};
};

Stats dpcp_decide_feasibility_ilp(DPCPInst& inputDpcp, Col& col, int timeLimit,
                                  std::ostream& log) {
  // Initial time instant
  auto start = std::chrono::high_resolution_clock::now();

  // Remove edges whose endpoints do not belongs to the same Pi and Qj
  Graph::edge_iterator ei, ei_end, next;
  boost::tie(ei, ei_end) = edges(g);
  for (next = ei; ei != ei_end; ei = next) {
    ++next;
    Vertex u = source(*ei, g);
    Vertex v = target(*ei, g);
    if (dpcp.get_P_part(u) != dpcp.get_P_part(v) &&
        dpcp.get_Q_part(u) != dpcp.get_Q_part(v)) {
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
    snprintf(name, sizeof(name), "x_%ld_%ld", dpcp.get_P_part(v),
             dpcp.get_Q_part(v));
    x[dpcp.get_current_id(v)] = IloBoolVar(cxenv, name);
  }

  // Objective
  IloExpr obj(cxenv);
  for (auto v : boost::make_iterator_range(vertices(g)))
    obj += x[dpcp.get_current_id(v)] * 1.0;
  cxobj = IloMaximize(cxenv, obj);
  cxmodel.add(cxobj);
  obj.end();

  // Constraints
  for (auto e : boost::make_iterator_range(edges(g))) {
    IloExpr expr(cxenv);
    expr += x[dpcp.get_current_id(source(e, g))] +
            x[dpcp.get_current_id(target(e, g))];
    cxcons.add(expr <= 1);
    expr.end();
  }
  cxmodel.add(cxcons);

  // Set parameters and extract the model
  cplex.extract(cxmodel);
  cplex.setOut(log);
  cplex.setDefaults();
  cplex.setParam(IloCplex::Param::Parallel, 1);  // Deterministic mode
  cplex.setParam(IloCplex::Param::Threads, 1);   // Single thread
  cplex.setParam(IloCplex::Param::TimeLimit, timeLimit);

  // Create the callback object
  EarlyStopCallback cb(dpcp.get_nP(), log);
  // Now we get to setting up the generic callback.
  // We instantiate a ThresholdCallback and set the contextMask parameter.
  CPXLONG contextMask = IloCplex::Callback::Context::Id::GlobalProgress;
  // If contextMask is not zero we add the callback.
  if (contextMask != 0) cplex.use(&cb, contextMask);

  // Finally, we can solve the model
  cplex.solve();

  // If optimal and obj value >= nP, then we have a feasible solution of DPCP
  // We must find the coloring in the original graph
  if (cplex.getCplexStatus() == IloCplex::Optimal &&
      static_cast<size_t>(cplex.getObjValue()) >= dpcp.get_nP()) {
    // Get variable values
    IloNumArray vals(cxenv);
    cplex.getValues(vals, x);
    // First, find selected vertices
    VertexVector selected;
    std::map<size_t, std::set<size_t>> adj;
    for (auto vv : boost::make_iterator_range(vertices(g))) {
      Vertex v = vertex(dpcp.get_current_id(vv), inputDpcp.get_graph());
      if (vals[dpcp.get_current_id(v)] > 0.5) {
        // Add adjacencies
        for (auto u : selected) {
          if (edge(u, v, inputDpcp.get_graph()).second) {
            size_t qj1 = inputDpcp.get_Q_part(u);
            size_t qj2 = inputDpcp.get_Q_part(v);
            if (qj1 < qj2)
              adj[qj1].insert(qj2);
            else if (qj2 > qj1)
              adj[qj2].insert(qj1);
          }
        }
        // Add new selected vertex
        selected.push_back(v);
      }
    }
    // Then, color them
    dpcp_dsatur_heur(inputDpcp, selected, adj, col);
    assert(col.check_coloring(inputDpcp));
    vals.end();
  }

  // Save stats
  Stats stats;
  if (cplex.getCplexStatus() == IloCplex::Optimal) {
    stats.state = static_cast<size_t>(cplex.getObjValue()) >= dpcp.get_nP()
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
  for (size_t v = 0; v < num_vertices(g); ++v) x[v].end();
  cxcons.end();
  cxmodel.end();
  cplex.end();
  cxenv.end();

  return stats;
}
