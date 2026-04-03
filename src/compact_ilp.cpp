#include "compact_ilp.hpp"

#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>

#include <cfloat>
#include <iostream>
#include <map>
#include <string>

#include "heur.hpp"

Stats solve_ilp(DPCPInst& dpcp, const Params& params, std::ostream& log,
                std::ostream& debugLog, Col& col) {
  Stats stats;

  // Try to find an initial coloring with the heuristic
  Col initialCol;
  if (params.heuristicRootNode == 1)
    stats = dpcp_1_step_greedy_heur(dpcp, initialCol);
  else if (params.heuristicRootNode == 2)
    stats = dpcp_2_step_greedy_heur(dpcp, initialCol, params);
  else if (params.heuristicRootNode == 3)
    stats = dpcp_2_step_semigreedy_heur(dpcp, initialCol, params);

  // Save initial solution stats
  if (params.heuristicRootNode >= 1 && params.heuristicRootNode <= 4) {
    stats.rootub = initialCol.get_n_colors();
    stats.rootHeurTime = stats.time;
  }

  // Number of colors
  size_t ncolors;
  if (initialCol.get_n_colors() > 0) {
    ncolors = initialCol.get_n_colors();
    log << "Initial coloring with " << initialCol.get_n_colors()
        << " colors found by heuristic." << std::endl;
  } else {
    ncolors = std::min(dpcp.get_nP(), dpcp.get_nQ());
    log << "No initial coloring found by heuristic." << std::endl;
  }

  // Initialize cplex enviroment
  IloEnv cxenv;
  IloModel cxmodel(cxenv);
  IloArray<IloNumVarArray> x(cxenv, num_vertices(dpcp.get_graph()));
  IloNumVarArray w(cxenv, ncolors);
  IloConstraintArray cxcons(cxenv);

  // Define variables
  for (size_t v = 0; v < num_vertices(dpcp.get_graph()); ++v) {
    x[v] = IloNumVarArray(cxenv, ncolors);
    for (size_t k = 0; k < ncolors; ++k) {
      char name[100];
      snprintf(name, sizeof(name), "x_%ld_%ld", v, k);
      x[v][k] = IloBoolVar(cxenv, name);
    }
  }
  for (size_t k = 0; k < ncolors; ++k) {
    char name[100];
    snprintf(name, sizeof(name), "w_%ld", k);
    w[k] = IloBoolVar(cxenv, name);
  }

  // Define objective
  IloExpr fobj(cxenv, 0);
  for (size_t k = 0; k < ncolors; ++k) fobj += w[k];
  cxmodel.add(IloMinimize(cxenv, fobj));

  // Constraints

  // \sum_{(a,b) \in V} \sum_{k \in C} x_(a,b)_k \geq 1, forall a \in A
  for (size_t pi = 0; pi < dpcp.get_nP(); ++pi) {
    auto& vec = dpcp.get_P()[pi];
    IloExpr restr(cxenv);
    for (Vertex v : vec)
      for (size_t k = 0; k < ncolors; ++k)
        restr += x[dpcp.get_current_id(v)][k];
    cxcons.add(restr >= 1);
  }

  // x_(a,b)_k + x_(a',b)_k' \leq 1, forall (a,b),(a',b) \in V with a != a',
  //                                         k, k' \in C with k != k'
  for (auto v1 : boost::make_iterator_range(vertices(dpcp.get_graph()))) {
    size_t pi1 = dpcp.get_P_part(v1);
    size_t qj1 = dpcp.get_Q_part(v1);
    size_t id1 = dpcp.get_current_id(v1);
    for (auto v2 : boost::make_iterator_range(vertices(dpcp.get_graph()))) {
      size_t pi2 = dpcp.get_P_part(v2);
      size_t qj2 = dpcp.get_Q_part(v2);
      size_t id2 = dpcp.get_current_id(v2);
      if ((qj1 != qj2) || (pi1 == pi2)) continue;
      for (size_t k1 = 0; k1 < ncolors; ++k1)
        for (size_t k2 = 0; k2 < ncolors; ++k2) {
          if (k1 == k2) continue;
          IloExpr restr(cxenv);
          restr += x[id1][k1] + x[id2][k2];
          cxcons.add(restr <= 1);
        }
    }
  }

  // x_(a,b)_k + x_(a',b')_k \leq w_k, forall ((a,b),(a',b')) \in E, k \in C
  for (auto e : boost::make_iterator_range(edges(dpcp.get_graph()))) {
    auto u = source(e, dpcp.get_graph());
    auto v = target(e, dpcp.get_graph());
    for (size_t k = 0; k < ncolors; ++k) {
      IloExpr restr(cxenv);
      restr +=
          x[dpcp.get_current_id(u)][k] + x[dpcp.get_current_id(v)][k] - w[k];
      cxcons.add(restr <= 0);
    }
  }

  cxmodel.add(cxcons);

  // Solve model
  IloCplex cplex(cxmodel);

  if (initialCol.get_n_colors() > 0) {
    // Mipstart
    Coloring coloring = initialCol.get_coloring();
    ColorClass classes = initialCol.get_color_classes();
    IloNumVarArray startVar(cxenv);
    IloNumArray startVal(cxenv);
    for (auto [idv, k] : coloring) {
      startVar.add(x[idv][k]);
      startVal.add(1);
    }
    for (auto& [k, s] : classes) {
      startVar.add(w[k]);
      startVal.add(1);
    }
    cplex.addMIPStart(startVar, startVal);
    startVal.end();
    startVar.end();
  }

  // Set parameters
  cplex.setDefaults();
  cplex.setOut(debugLog);
  cplex.setParam(IloCplex::Param::TimeLimit, params.timeLimit);
  cplex.setParam(IloCplex::Param::Parallel, 1);  // Deterministic mode
  cplex.setParam(IloCplex::Param::Threads, 1);   // Single thread
  // cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicEffort, 0);

  // Solve
  cplex.solve();

  // Get final state
  STATE state;
  switch (cplex.getCplexStatus()) {
    case IloCplex::CplexStatus::Optimal:
      state = OPTIMAL;
      break;
    case IloCplex::CplexStatus::Infeasible:
      state = INFEASIBLE;
      break;
    case IloCplex::CplexStatus::AbortTimeLim:
      if (cplex.getSolnPoolNsolns())
        state = FEASIBLE;
      else
        state = TIME_EXCEEDED;
      break;
    case IloCplex::CplexStatus::MemLimFeas:
    case IloCplex::CplexStatus::MemLimInfeas:
      state = MEM_EXCEEDED;
      break;
    default:
      state = UNKNOWN;
      break;
  }

  if (state == OPTIMAL || state == FEASIBLE) {
    // Recover coloring
    for (auto v : boost::make_iterator_range(vertices(dpcp.get_graph())))
      for (size_t k = 0; k < ncolors; ++k)
        if (cplex.getValue(x[dpcp.get_current_id(v)][k]) > 0.5)
          col.set_color(dpcp, dpcp.get_current_id(v), k);
    assert(col.check_coloring(dpcp));
  }

  // Complete stats
  stats.nvars = static_cast<int>(cplex.getNcols());
  stats.ncons = static_cast<int>(cplex.getNrows());
  stats.state = state;
  stats.time = cplex.getTime();
  stats.nodes = static_cast<int>(cplex.getNnodes());
  stats.lb = cplex.getBestObjValue();
  stats.ub = -1;
  if (state == OPTIMAL || state == FEASIBLE) {
    stats.ub = static_cast<int>(cplex.getObjValue() + 0.5);
    stats.gap = cplex.getMIPRelativeGap();
  }

  // Free memory
  fobj.end();
  cxcons.end();
  for (size_t v = 0; v < num_vertices(dpcp.get_graph()); ++v) x[v].end();
  x.end();
  w.end();
  cplex.end();
  cxmodel.end();
  cxenv.end();

  return stats;
}
