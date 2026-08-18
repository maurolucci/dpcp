#include "compact_ilp.hpp"

#include <gurobi_c++.h>
#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>

#include <cfloat>
#include <chrono>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "heur.hpp"

Stats solve_ilp_cplex(DPCPInst& _dpcp, const Params& params, std::ostream& log,
                      std::ostream& debugLog, Col& _col) {
  Stats stats;
  HeurStats heurStats;

  // Initial time
  auto startTime = std::chrono::high_resolution_clock::now();

  // Make a copy of dpcp to avoid modifying the original instance
  DPCPInst dpcp(_dpcp);

  // Preprocess the instance
  dpcp.preprocess(true, params.preprocessing);
  if (dpcp.is_infeasible_instance()) {
    stats.state = INFEASIBLE;
    stats.time = std::chrono::duration<double>(
                     std::chrono::high_resolution_clock::now() - startTime)
                     .count();
    return stats;
  }

  // Try to find an initial coloring with the heuristic
  Col initialCol;
  if (params.heuristicInitial == 1)
    heurStats = dpcp_1_step_greedy_heur(dpcp, initialCol);
  else if (params.heuristicInitial == 2)
    heurStats = dpcp_2_step_greedy_heur(dpcp, initialCol, params);
  else if (params.heuristicInitial == 3)
    heurStats = dpcp_2_step_semigreedy_heur(dpcp, initialCol, params);

  // Save initial solution stats
  if (params.heuristicInitial >= 1 && params.heuristicInitial <= 4) {
    stats.initialHeurValue = initialCol.get_n_colors();
    stats.initialHeurTime = heurStats.totalTime;
    stats.initialSemigreedyIters = (params.heuristicInitial == 3)
                                       ? static_cast<int>(heurStats.totalIters)
                                       : 0;
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
  cplex.setParam(IloCplex::Param::TimeLimit,
                 params.timeLimit -
                     std::chrono::duration<double>(
                         std::chrono::high_resolution_clock::now() - startTime)
                         .count());
  cplex.setParam(IloCplex::Param::Parallel, 1);      // Deterministic mode
  cplex.setParam(IloCplex::Param::Threads, 1);       // Single thread
  cplex.setParam(IloCplex::Param::WorkMem, 4096.0);  // Memory limit in MB
  cplex.setParam(IloCplex::Param:: ::Limits::TreeMemory,
                 5120.0);  // Tree memory limit in MB
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

  if (state == OPTIMAL || ((state == TIME_EXCEEDED || state == MEM_EXCEEDED) &&
                           cplex.getSolnPoolNsolns())) {
    // Recover coloring
    Col col;
    for (auto v : boost::make_iterator_range(vertices(dpcp.get_graph())))
      for (size_t k = 0; k < ncolors; ++k)
        if (cplex.getValue(x[dpcp.get_current_id(v)][k]) > 0.5)
          col.set_color(dpcp, dpcp.get_current_id(v), k);
    assert(col.check_coloring(dpcp));

    // Translate coloring to original instance
    _col = col.translate_coloring(dpcp, _dpcp);
    assert(_col.check_coloring(_dpcp));
  }

  // Complete stats
  stats.nvars = static_cast<int>(cplex.getNcols());
  stats.ncons = static_cast<int>(cplex.getNrows());
  stats.state = state;
  stats.time = std::chrono::duration<double>(
                   std::chrono::high_resolution_clock::now() - startTime)
                   .count();
  stats.nodes = static_cast<int>(cplex.getNnodes());
  stats.lb = cplex.getBestObjValue();
  stats.ub = -1;
  if (state == OPTIMAL || ((state == TIME_EXCEEDED || state == MEM_EXCEEDED) &&
                           cplex.getSolnPoolNsolns())) {
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

class StreamLogger : public GRBCallback {
 private:
  std::ostream& stream_;

 protected:
  void callback() override {
    try {
      if (where == GRB_CB_MESSAGE) {
        // Retrieve the log message string from Gurobi
        std::string msg = getStringInfo(GRB_CB_MSG);
        stream_ << msg;  // Stream to your ostream destination
        stream_.flush();
      }
    } catch (...) {
      // Suppress or handle exceptions in callback safely
    }
  }

 public:
  StreamLogger(std::ostream& out) : stream_(out) {}
};

Stats solve_ilp_gurobi(DPCPInst& _dpcp, const Params& params, std::ostream& log,
                       std::ostream& debugLog, Col& _col) {
  Stats stats;
  HeurStats heurStats;
  auto startTime = std::chrono::high_resolution_clock::now();
  DPCPInst dpcp(_dpcp);

  dpcp.preprocess(true, params.preprocessing);
  if (dpcp.is_infeasible_instance()) {
    stats.state = INFEASIBLE;
    stats.time = std::chrono::duration<double>(
                     std::chrono::high_resolution_clock::now() - startTime)
                     .count();
    return stats;
  }

  Col initialCol;
  if (params.heuristicInitial == 1)
    heurStats = dpcp_1_step_greedy_heur(dpcp, initialCol);
  else if (params.heuristicInitial == 2)
    heurStats = dpcp_2_step_greedy_heur(dpcp, initialCol, params);
  else if (params.heuristicInitial == 3)
    heurStats = dpcp_2_step_semigreedy_heur(dpcp, initialCol, params);

  if (params.heuristicInitial >= 1 && params.heuristicInitial <= 4) {
    stats.initialHeurValue = initialCol.get_n_colors();
    stats.initialHeurTime = heurStats.totalTime;
    stats.initialSemigreedyIters = (params.heuristicInitial == 3)
                                       ? static_cast<int>(heurStats.totalIters)
                                       : 0;
  }

  size_t ncolors;
  if (initialCol.get_n_colors() > 0) {
    ncolors = initialCol.get_n_colors();
    log << "Initial coloring with " << initialCol.get_n_colors()
        << " colors found by heuristic." << std::endl;
  } else {
    ncolors = std::min(dpcp.get_nP(), dpcp.get_nQ());
    log << "No initial coloring found by heuristic." << std::endl;
  }

  GRBEnv env(true);
  env.start();
  GRBModel model(env);
  std::vector<std::vector<GRBVar>> x(num_vertices(dpcp.get_graph()),
                                     std::vector<GRBVar>(ncolors));
  std::vector<GRBVar> w(ncolors);

  for (size_t v = 0; v < num_vertices(dpcp.get_graph()); ++v)
    for (size_t k = 0; k < ncolors; ++k)
      x[v][k] =
          model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                       "x_" + std::to_string(v) + "_" + std::to_string(k));
  for (size_t k = 0; k < ncolors; ++k)
    w[k] = model.addVar(0.0, 1.0, 1.0, GRB_BINARY, "w_" + std::to_string(k));

  GRBLinExpr objective = 0;
  for (size_t k = 0; k < ncolors; ++k) objective += w[k];
  model.setObjective(objective, GRB_MINIMIZE);

  for (size_t pi = 0; pi < dpcp.get_nP(); ++pi) {
    GRBLinExpr restr = 0;
    for (Vertex v : dpcp.get_P()[pi])
      for (size_t k = 0; k < ncolors; ++k)
        restr += x[dpcp.get_current_id(v)][k];
    model.addConstr(restr >= 1);
  }

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
        for (size_t k2 = 0; k2 < ncolors; ++k2)
          if (k1 != k2) model.addConstr(x[id1][k1] + x[id2][k2] <= 1);
    }
  }

  for (auto e : boost::make_iterator_range(edges(dpcp.get_graph()))) {
    auto u = source(e, dpcp.get_graph());
    auto v = target(e, dpcp.get_graph());
    for (size_t k = 0; k < ncolors; ++k)
      model.addConstr(
          x[dpcp.get_current_id(u)][k] + x[dpcp.get_current_id(v)][k] <= w[k]);
  }

  if (initialCol.get_n_colors() > 0) {
    for (auto [idv, k] : initialCol.get_coloring())
      x[idv][k].set(GRB_DoubleAttr_Start, 1.0);
    for (auto& [k, _] : initialCol.get_color_classes())
      w[k].set(GRB_DoubleAttr_Start, 1.0);
  }

  // Instantiate your ostream destination for logging
  model.set(GRB_IntParam_LogToConsole, 0);
  StreamLogger my_logger(debugLog);
  model.setCallback(&my_logger);
  const double elapsed =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now() -
                                    startTime)
          .count();
  model.set(GRB_DoubleParam_TimeLimit,
            std::max(0.0, static_cast<double>(params.timeLimit) - elapsed));
  model.set(GRB_IntParam_Threads, 1);
  model.set(GRB_DoubleParam_NodefileStart, 4.0);
  model.set(GRB_DoubleParam_SoftMemLimit, 5.0);
  model.optimize();

  const int status = model.get(GRB_IntAttr_Status);
  STATE state = UNKNOWN;
  if (status == GRB_OPTIMAL)
    state = OPTIMAL;
  else if (status == GRB_INFEASIBLE)
    state = INFEASIBLE;
  else if (status == GRB_TIME_LIMIT)
    state = TIME_EXCEEDED;
  else if (status == GRB_MEM_LIMIT)
    state = MEM_EXCEEDED;

  const bool hasSolution = model.get(GRB_IntAttr_SolCount) > 0;
  if (state == OPTIMAL ||
      ((state == TIME_EXCEEDED || state == MEM_EXCEEDED) && hasSolution)) {
    Col col;
    for (auto v : boost::make_iterator_range(vertices(dpcp.get_graph())))
      for (size_t k = 0; k < ncolors; ++k)
        if (x[dpcp.get_current_id(v)][k].get(GRB_DoubleAttr_X) > 0.5)
          col.set_color(dpcp, dpcp.get_current_id(v), k);
    assert(col.check_coloring(dpcp));
    _col = col.translate_coloring(dpcp, _dpcp);
    assert(_col.check_coloring(_dpcp));
  }

  stats.nvars = model.get(GRB_IntAttr_NumVars);
  stats.ncons = model.get(GRB_IntAttr_NumConstrs);
  stats.state = state;
  stats.time = std::chrono::duration<double>(
                   std::chrono::high_resolution_clock::now() - startTime)
                   .count();
  stats.nodes = static_cast<int>(model.get(GRB_DoubleAttr_NodeCount));
  stats.lb = model.get(GRB_DoubleAttr_ObjBound);
  stats.ub = -1;
  if (state == OPTIMAL ||
      ((state == TIME_EXCEEDED || state == MEM_EXCEEDED) && hasSolution)) {
    stats.ub = static_cast<int>(model.get(GRB_DoubleAttr_ObjVal) + 0.5);
    stats.gap = model.get(GRB_DoubleAttr_MIPGap);
  }
  debugLog.flush();
  return stats;
}
