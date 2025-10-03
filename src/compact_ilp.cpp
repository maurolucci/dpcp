#include "compact_ilp.hpp"
#include "heur.hpp"

#include <cfloat>
#include <iostream>
#include <map>
#include <string>

#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>

Stats solve_ilp(const Graph &graph, const Params &params, std::ostream &log,
                Col &col) {

  Stats stats;

  // Try to find an initial coloring with the heuristic
  Col initialCol;
  GraphEnv genv(graph, params, true);
  if (params.heuristicRootNode == 1)
    stats = dpcp_1_step_greedy_heur(genv, initialCol);
  else if (params.heuristicRootNode == 2)
    stats =
        dpcp_1_step_semigreedy_heur(genv, initialCol, params.heuristicRootIter);
  else if (params.heuristicRootNode == 3)
    stats = dpcp_2_step_greedy_heur(genv, initialCol);
  else if (params.heuristicRootNode == 4)
    stats =
        dpcp_2_step_semigreedy_heur(genv, initialCol, params.heuristicRootIter);

  // Save initial solution stats
  if (params.heuristicRootNode >= 1 && params.heuristicRootNode <= 4) {
    stats.initSol = initialCol.get_n_colors();
    stats.initSolTime = stats.time;
  }

  // Number of colors
  size_t ncolors;
  if (initialCol.get_n_colors() > 0) {
    ncolors = initialCol.get_n_colors();
    log << "Initial coloring with " << initialCol.get_n_colors()
        << " colors found by heuristic." << std::endl;
  } else {
    ncolors = std::min(genv.nA, genv.nB);
    log << "No initial coloring found by heuristic." << std::endl;
  }

  // Initialize cplex enviroment
  IloEnv cxenv;
  IloModel cxmodel(cxenv);
  IloArray<IloNumVarArray> x(cxenv, num_vertices(graph));
  IloNumVarArray w(cxenv, ncolors);
  IloConstraintArray cxcons(cxenv);

  // Define variables
  for (size_t v = 0; v < num_vertices(graph); ++v) {
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
  for (size_t k = 0; k < ncolors; ++k)
    fobj += w[k];
  cxmodel.add(IloMinimize(cxenv, fobj));

  // Constraints

  // \sum_{(a,b) \in V} \sum_{k \in C} x_(a,b)_k \geq 1, forall a \in A
  for (auto &[a, vec] : genv.Va) {
    IloExpr restr(cxenv);
    for (Vertex v : vec)
      for (size_t k = 0; k < ncolors; ++k)
        restr += x[graph[v].id][k];
    cxcons.add(restr >= 1);
  }

  // x_(a,b)_k + x_(a',b)_k' \leq 1, forall (a,b),(a',b) \in V with a != a',
  //                                         k, k' \in C with k != k'
  for (auto v1 : boost::make_iterator_range(vertices(graph))) {
    auto [a1, b1, id1] = graph[v1];
    for (auto v2 : boost::make_iterator_range(vertices(graph))) {
      auto [a2, b2, id2] = graph[v2];
      if ((b1 != b2) || (a1 == a2))
        continue;
      for (size_t k1 = 0; k1 < ncolors; ++k1)
        for (size_t k2 = 0; k2 < ncolors; ++k2) {
          if (k1 == k2)
            continue;
          IloExpr restr(cxenv);
          restr += x[id1][k1] + x[id2][k2];
          cxcons.add(restr <= 1);
        }
    }
  }

  // x_(a,b)_k + x_(a',b')_k \leq w_k, forall ((a,b),(a',b')) \in E, k \in C
  for (auto e : boost::make_iterator_range(edges(graph))) {
    auto u = source(e, graph);
    auto v = target(e, graph);
    for (size_t k = 0; k < ncolors; ++k) {
      IloExpr restr(cxenv);
      restr += x[graph[u].id][k] + x[graph[v].id][k] - w[k];
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
    for (auto [v, k] : coloring) {
      startVar.add(x[graph[v].id][k]);
      startVal.add(1);
    }
    for (auto &[k, s] : classes) {
      startVar.add(w[k]);
      startVal.add(1);
    }
    cplex.addMIPStart(startVar, startVal);
    startVal.end();
    startVar.end();
    stats.initSol = initialCol.get_n_colors();
  }

  // Set parameters
  cplex.setDefaults();
  cplex.setOut(log);
  cplex.setParam(IloCplex::Param::TimeLimit, params.timeLimit);
  cplex.setParam(IloCplex::Param::Parallel, 1); // Deterministic mode
  cplex.setParam(IloCplex::Param::Threads, 1);  // Single thread
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
    for (auto v : boost::make_iterator_range(vertices(graph)))
      for (size_t k = 0; k < ncolors; ++k)
        if (cplex.getValue(x[graph[v].id][k]) > 0.5)
          col.set_color(graph, v, k);
    assert(col.check_coloring(graph));
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

  return stats;
}