#include "lp.hpp"

#include <cfloat>
#include <limits>
#include <numeric>
#include <queue>

#include "feas.hpp"
#include "heur.hpp"
#include "pricing.hpp"
#include "random.hpp"

double get_elapsed_time(auto startTime) {
  return std::chrono::duration_cast<std::chrono::duration<double>>(
             std::chrono::high_resolution_clock::now() - startTime)
      .count();
}

auto find_most_fractional(std::map<Vertex, double>& m) {
  return std::max_element(m.begin(), m.end(),
                          [](const std::pair<Vertex, double>& a,
                             const std::pair<Vertex, double>& b) -> bool {
                            return std::abs(a.second - 0.5) >
                                   std::abs(b.second - 0.5);
                          });
}

LP::LP(DPCPInst&& dpcp, Pool&& pool, const DPCPInst& origDpcp, Params& params,
       Stats& stats, std::ostream& log, std::ostream& debugLog, bool isRoot)
    : dpcp(std::move(dpcp)),
      pool(std::move(pool)),
      origDpcp(origDpcp),
      params(params),
      stats(stats),
      log(log),
      debugLog(debugLog),
      isRoot(isRoot),
      objVal(-1.0),
      state(LP_UNSOLVED),
      integerSource(LP_INTEGER_SOURCE_NONE),
      initializedWithDummy(false),
      stables(),
      posVars() {}

LP::LP(const LP& other)
    : dpcp(other.dpcp),
      pool(),
      lateColumns(),
      origDpcp(other.origDpcp),
      params(other.params),
      stats(other.stats),
      log(other.log),
      debugLog(other.debugLog),
      isRoot(false),
      objVal(-1.0),
      state(LP_UNSOLVED),
      integerSource(LP_INTEGER_SOURCE_NONE),
      initializedWithDummy(false),
      stables(),
      posVars() {}

LP::LP(LP&& other) noexcept
    : dpcp(std::move(other.dpcp)),
      pool(std::move(other.pool)),
      lateColumns(std::move(other.lateColumns)),
      origDpcp(other.origDpcp),
      params(other.params),
      stats(other.stats),
      log(other.log),
      debugLog(other.debugLog),
      isRoot(other.isRoot),
      objVal(other.objVal),
      state(other.state),
      integerSource(other.integerSource),
      initializedWithDummy(other.initializedWithDummy),
      stables(std::move(other.stables)),
      posVars(std::move(other.posVars)),
      branchingVertex(other.branchingVertex),
      coloring(std::move(other.coloring)),
      pricingSummary(std::move(other.pricingSummary)) {}

LP::~LP() {}

LP_STATE LP::solve(double timelimit, double ub) {
  auto startTime = std::chrono::high_resolution_clock::now();

  // Cleaning
  initializedWithDummy = false;
  stables.clear();
  posVars.clear();
  coloring.reset_coloring();
  state = LP_UNSOLVED;
  integerSource = LP_INTEGER_SOURCE_NONE;

  // Infeasibility check
  if (dpcp.is_infeasible_instance()) {
    stats.ninfeasPrepro++;
    state = LP_INFEASIBLE;
    if (params.is_verbose(2))
      debugLog << "LP detected infeasibility." << std::endl;
    return state;
  }

  // Trivial solution check
  if (dpcp.has_trivial_solution()) {
    stats.ntrivial++;
    objVal = 1.0;
    state = LP_INTEGER;
    integerSource = LP_INTEGER_SOURCE_TRIVIAL;
    if (params.is_verbose(2))
      debugLog << "LP found trivial solution." << std::endl;
    return state;
  }

  // GCP instance check
  if (dpcp.is_gcp_instance()) {
    if (params.is_verbose(2))
      debugLog << "LP reduced to GCP instance." << std::endl;
    return gcp_solve(timelimit, ub);
  }

  // Apply heuristic at the current node
  heuristic_solve();
  if (params.is_verbose(2)) {
    if (has_heur_solution())
      debugLog << "LP heuristic solution: " << coloring.get_n_colors()
               << " colors, " << get_elapsed_time(startTime) << " seconds."
               << std::endl;
    else
      debugLog << "LP heuristic failed to find a solution." << std::endl;
  }

  // Apply feasibility check at the current node
  if (!feasibility_solve()) {
    stats.ninfeasCheck++;
    state = LP_INFEASIBLE;
    if (params.is_verbose(2))
      debugLog << "LP feasibility check proved infeasibility." << std::endl;
    return state;
  }

  // Initialize cplex environment
  CplexEnv cenv;
  IloCplex cplex(cenv.Xmodel);
  set_parameters(cenv, cplex);
  add_constraints_and_objective(cenv);

  // Add initial columns
  add_initial_columns(cenv);

  // Add late columns (e.g., positive parent columns in mode 3)
  for (auto& col : lateColumns) {
    add_column(cenv, col);
  }

  if (params.is_verbose(2)) {
    debugLog << "LP initialized with " << stables.size() << " columns ("
             << lateColumns.size() << " late), using "
             << (initializedWithDummy ? "dummy column" : "initial solution")
             << "." << std::endl;
  }

  // Initialize arrays for dual values
  IloNumArray dualsP(cenv.Xenv, dpcp.get_nP());
  IloNumArray dualsQ(cenv.Xenv, dpcp.get_nQ());

  // Initialize pricing environment
  PricingEnv penv(dpcp, params.pricingExactTimeLimit);

  // Reset per-LP pricing summary
  pricingSummary = PricingSummary();
  pricingSummary.iters = 0;

  // Start the column generation loop
  while (state == LP_UNSOLVED) {
    // Set time limit
    double timelimit2 = timelimit - get_elapsed_time(startTime);
    if (timelimit2 < 0) {
      state = LP_TIME_EXCEEDED;
      break;
    }
    cplex.setParam(IloCplex::Param::TimeLimit, timelimit2);

    // CPLEX optimization
    cplex.solve();

    // Handle errors
    IloCplex::CplexStatus status = cplex.getCplexStatus();
    if (status == IloCplex::AbortTimeLim) {
      state = LP_TIME_EXCEEDED;
      break;
    } else if (status == IloCplex::MemLimFeas ||
               status == IloCplex::MemLimInfeas) {
      state = LP_MEM_EXCEEDED;
      break;
    }

    // Compute dual value array
    cplex.getDuals(dualsP, cenv.XrestrP);
    cplex.getDuals(dualsQ, cenv.XrestrQ);
    objVal = cplex.getObjValue();

    // Pricing
    pricingSummary.iters++;
    int ret = pricing(cenv, penv, dualsP, dualsQ);

    if (ret > 0)
      continue;  // Keep generating columns
    else if (ret == 0)
      break;  // Optimality reached
    else if (ret == -1) {
      state = LP_TIME_EXCEEDED_PR;
      break;
    } else if (ret == -2) {
      state = LP_MEM_EXCEEDED_PR;
      break;
    } else {
      state = LP_UNKNOWN;
      break;
    }
  }

  // Update global stats once using the accumulated pricing summary.
  update_stats_from_pricing_summary();
  log_pricing_summary();

  // Save root lower bound
  if (isRoot) {
    stats.rootlb = objVal;
  }

  // If no error occurred, recover primal values and objective value to check
  // integrality
  if (state == LP_UNSOLVED) {
    IloNumArray values = IloNumArray(cenv.Xenv, cenv.Xvars.getSize());
    cplex.getValues(values, cenv.Xvars);

    // Check if dummy column was used
    if (initializedWithDummy && values[0] > EPSILON) {
      // If the objective value is greater than W, then the instance is
      // infeasible
      if (objVal > std::min(dpcp.get_nP(), dpcp.get_nQ()) + EPSILON) {
        state = LP_INFEASIBLE;
        stats.ninfeasAux++;
      }
      // Otherwise, the initialization failed
      else {
        log << "Initialization failed" << std::endl;
        log << "Dummy variable has value " << values[0]
            << " and the objective value is " << objVal
            << " <= " << std::min(dpcp.get_nP(), dpcp.get_nQ()) + EPSILON
            << std::endl;
        state = LP_INIT_FAIL;
      }
    }

    // Initial variable
    int i = initializedWithDummy ? 1 : 0;

    // Integrality check
    if (state == LP_UNSOLVED) {
      state = LP_INTEGER;
      for (; i < cenv.Xvars.getSize(); ++i) {
        if (values[i] < EPSILON)
          continue;
        else if (values[i] < 1 - EPSILON)
          state = LP_FRACTIONAL;
        posVars.push_back(i);
      }

      if (state == LP_FRACTIONAL) {
        // Find branching variable
        branchingVertex = params.use_fms_branching()
                              ? get_branching_variable_FMS(values)
                              : get_branching_variable_LNTT(values);
      } else if (state == LP_INTEGER) {
        objVal = posVars.size();
        integerSource = LP_INTEGER_SOURCE_LR;
      }
    }

    if (params.is_verbose(2)) {
      debugLog << "LP solve end: state=" << state << "("
               << get_lp_state_as_str(state) << "), objVal=" << objVal
               << ", #posVars=" << posVars.size()
               << ", time=" << get_elapsed_time(startTime) << std::endl;
    }

    values.end();
  }

  cplex.end();
  return state;
}

// Solve a graph coloring problem instance with exactcolors
LP_STATE LP::gcp_solve(double timelimit, double ub) {
  if (timelimit <= 0) {
    state = LP_TIME_EXCEEDED;
    return state;
  }

  auto startTime = std::chrono::high_resolution_clock::now();

  COLORproblem colorproblem;
  COLORparms* parms = &(colorproblem.parms);
  colordata* cd = &(colorproblem.root_cd);
  int ncolors = 0;
  const Partition& Q = dpcp.get_Q();
  size_t nP = dpcp.get_nP();
  size_t nQ = dpcp.get_nQ();

  // Build GCP graph
  GCPGraph graphGCP = dpcp.get_gcp_graph();

  // Get edge list
  int ecount = 0;
  int elist[2 * num_edges(graphGCP)];
  for (auto e : boost::make_iterator_range(edges(graphGCP))) {
    elist[2 * ecount] = source(e, graphGCP);
    elist[2 * ecount++ + 1] = target(e, graphGCP);
  }

  // Initialization
  COLORset* colorclasses = (COLORset*)NULL;
  COLORproblem_init_with_graph(&colorproblem, num_vertices(graphGCP), ecount,
                               elist);
  cd->id = 0;
  colorproblem.ncolordata = 1;
  parms->branching_cpu_limit = timelimit;
  colorproblem.root_cd.upper_bound = (ub < std::numeric_limits<double>::max())
                                         ? static_cast<int>(ub)
                                         : std::min(nP, nQ);  // UB
  COLORexact_coloring(&colorproblem, &ncolors, &colorclasses);

  // Optimality check
  if (cd->lower_bound == cd->upper_bound) {
    state = LP_INTEGER;
    integerSource = LP_INTEGER_SOURCE_GCP;
    objVal = ncolors;
    // Save solution
    // Stable sets need to be translated into vertices of the original graph
    for (int i = 0; i < ncolors; ++i) {
      colorclasses[i].age = 0;

      // Find the correct size
      int newCount = 0;
      for (int j = 0; j < colorclasses[i].count; ++j)
        newCount += Q[colorclasses[i].members[j]].size();

      Column stable;

      // Write translated stable set
      for (int j = 0; j < colorclasses[i].count; ++j)
        for (auto v : Q[colorclasses[i].members[j]])
          stable.add_vertex(v, dpcp.get_P_part(v), colorclasses[i].members[j]);

      // Add the translated stable set
      stables.push_back(stable);
      posVars.push_back(i);
    }
  } else
    state = LP_TIME_EXCEEDED;

  COLORproblem_free(&colorproblem);
  COLORfree_sets(&colorclasses, &ncolors);
  COLORlp_free_env();
  stats.gcpTime += std::chrono::duration<double>(
                       std::chrono::high_resolution_clock::now() - startTime)
                       .count();
  return state;
}

// Heuristic solution of the DPCP instances at the current node
// The heuristic used is selected depending on the parameters
void LP::heuristic_solve() {
  HeurStats heurStats;
  int heur = isRoot ? params.heuristicRootNode : params.heuristicOtherNodes;
  switch (heur) {
    case 0:
      return;
    case 1:
      heurStats = dpcp_1_step_greedy_heur(dpcp, coloring);
      break;
    case 2:
      heurStats = dpcp_2_step_greedy_heur(dpcp, coloring, params);
      break;
    case 3:
      assert(dpcp.check_consistency());
      heurStats = dpcp_2_step_semigreedy_heur(dpcp, coloring, params);
      break;
    default:
      log << "Warning: unknown heuristic code " << heur
          << ", skipping heuristic solve." << std::endl;
      return;
  }

  // Update heuristic stats
  if (isRoot) {
    stats.rootHeurTime = heurStats.totalTime;
    if (heurStats.state == FEASIBLE) stats.rootub = heurStats.value;
  } else {
    stats.otherNodesHeurTime += heurStats.totalTime;
  }

  return;
}

// Feasibility check of the DPCP instance at the current node
// The feasibility check used is selected depending on the parameters
bool LP::feasibility_solve() {
  if (has_heur_solution()) return true;
  Stats feasStats;
  int check =
      isRoot ? params.feasibilityRootNode : params.feasibilityOtherNodes;
  if (check == 0) return true;
  int timeLimit = isRoot ? params.feasibilityRootNodeTimeLimit
                         : params.feasibilityOtherNodesTimeLimit;
  struct NullBuffer : std::streambuf {
    int overflow(int c) override { return c; }
  } nullBuffer;
  std::ostream nullstream(&nullBuffer);

  // Make a copy of the DPCP instance to avoid modifying the one at the current
  // node
  DPCPInst dpcpCopy(dpcp);
  if (check == 1)
    feasStats =
        dpcp_decide_feasibility_enumerative(dpcpCopy, coloring, nullstream);
  else if (check == 2)
    feasStats =
        dpcp_decide_feasibility_ilp(dpcpCopy, coloring, timeLimit, nullstream);
  else {
    log << "Warning: unknown feasibility check code " << check
        << ", skipping feasibility solve." << std::endl;
    return true;
  }

  // Update feasibility check stats
  if (isRoot)
    stats.rootFeasTime = feasStats.time;
  else {
    stats.otherNodesFeasNCalls++;
    stats.otherNodesFeasTime += feasStats.time;
  }

  return feasStats.state == FEASIBLE;
}

// Set CPLEX's parameters
void LP::set_parameters(CplexEnv& cenv, IloCplex& cplex) {
  cplex.setDefaults();
  cplex.setParam(IloCplex::Param::Threads, 1);
  cplex.setParam(IloCplex::Param::Parallel, 1);
  cplex.setOut(cenv.Xenv.getNullStream());
  cplex.setWarning(cenv.Xenv.getNullStream());
}

// Initialize the objective function and the constraints of the LP
void LP::add_constraints_and_objective(CplexEnv& cenv) {
  // Add ">= 1" constraints, one for each part in P
  for (auto pi = dpcp.get_nP(); pi > 0; --pi)
    cenv.XrestrP.add(IloRange(cenv.Xenv, 1.0, IloInfinity));
  // Add "<= 1" constraints, one for each part in Q
  for (auto qj = dpcp.get_nQ(); qj > 0; --qj)
    cenv.XrestrQ.add(IloRange(cenv.Xenv, -1.0, IloInfinity));
  cenv.Xmodel.add(cenv.XrestrP);
  cenv.Xmodel.add(cenv.XrestrQ);
  // Add objective function
  cenv.Xobj = IloMinimize(cenv.Xenv, 0.0);
  cenv.Xmodel.add(cenv.Xobj);
  return;
}

// Add initial columns to the LR from a heuristic solution if it exists, or a
// dummy column otherwise. The heuristic soluction may be found by the
// heuristic_solve() method or the feasibility_solve() method, depending on the
// parameters of the algorithm. The
void LP::add_initial_columns(CplexEnv& cenv) {
  // Try to initialize with the heuristic solution
  if (has_heur_solution()) {
    for (size_t k = 0; k < coloring.get_n_colors(); ++k) {
      Column stab = coloring.get_stable(dpcp, k);
      add_column(cenv, stab);
    }
  }
  // Otherwise, initialize with a dummy column
  else {
    // Add a new dummy variable z with z_i = 1 forall i in P and z_j = 0
    // forall j in Q and set a high cost for z
    IloNumColumn z = cenv.Xobj(params.initializationBigWeight);
    for (size_t pi = 0; pi < dpcp.get_nP(); ++pi) z += cenv.XrestrP[pi](1.0);
    for (size_t qj = 0; qj < dpcp.get_nQ(); ++qj) z += cenv.XrestrQ[qj](0.0);
    cenv.Xvars.add(IloNumVar(z));
    stables.push_back(Column());  // Push a null column
    initializedWithDummy = true;
  }
  return;
}

int LP::pricing(CplexEnv& cenv, PricingEnv& penv, IloNumArray& dualsP,
                IloNumArray& dualsQ) {
  int addedCols = 0;
  bool useGreedy = false;
  bool usePQmwss = false;
  bool usePmwss = false;
  int pricingOrder = 1;
  int pricingMethod = params.pricingMethod;

  // Resolve pricing method 6, which depends on the density of the graph
  if (pricingMethod == 6) {
    if (dpcp.get_density() <= 0.6) {
      pricingMethod = 2;
    } else {
      pricingMethod = 5;
    }
  }

  // Resolve pricing pipeline directly from pricingMethod.
  switch (pricingMethod) {
    case 0:
      break;
    case 1:
      useGreedy = true;
      break;
    case 2:
      useGreedy = true;
      usePQmwss = true;
      pricingOrder = 1;
      break;
    case 3:
      useGreedy = true;
      usePmwss = true;
      pricingOrder = 2;
      break;
    case 4:
      useGreedy = true;
      usePQmwss = true;
      usePmwss = true;
      pricingOrder = 1;
      break;
    case 5:
      useGreedy = true;
      usePQmwss = true;
      usePmwss = true;
      pricingOrder = 2;
      break;
    default:
      break;
  }

  // First, look for entering columns in the pool
  addedCols = pricing_pool(cenv, dualsP, dualsQ);
  if (addedCols > 0) return addedCols;

  // Second, look for entering columns with a greedy heuristic
  addedCols = pricing_greedy(cenv, penv, dualsP, dualsQ, useGreedy);
  if (addedCols > 0) return addedCols;

  if (pricingOrder == 1) {
    // Third, P,Q-MWSSP heuristic
    addedCols = pricing_P_Q_mwss(cenv, penv, dualsP, dualsQ, usePQmwss);
    if (addedCols > 0) return addedCols;

    // Fourth, P-MWSSP heuristic
    addedCols = pricing_P_mwss(cenv, penv, dualsP, dualsQ, usePmwss);
    if (addedCols >= 0) return addedCols;
  } else {
    // Third, P-MWSSP heuristic
    addedCols = pricing_P_mwss(cenv, penv, dualsP, dualsQ, usePmwss);
    if (addedCols >= 0) return addedCols;

    // Fourth, P,Q-MWSSP heuristic
    addedCols = pricing_P_Q_mwss(cenv, penv, dualsP, dualsQ, usePQmwss);
    if (addedCols > 0) return addedCols;
  }

  // Fifth, exact resolution of pricing
  return pricing_exact(cenv, penv, dualsP, dualsQ);
}

void LP::update_stats_from_pricing_summary() {
  if (isRoot) {
    stats.rootNCallsPool += static_cast<int>(pricingSummary.callsPool);
    stats.rootNCallsHeur += static_cast<int>(pricingSummary.callsGreedy);
    stats.rootNCallsMwis1 += static_cast<int>(pricingSummary.callsPQmwss);
    stats.rootNCallsMwis2 += static_cast<int>(pricingSummary.callsPmwss);
    stats.rootNCallsExact += static_cast<int>(pricingSummary.callsExact);

    stats.rootNColsPool += static_cast<int>(pricingSummary.colsPool);
    stats.rootNColsHeur += static_cast<int>(pricingSummary.colsGreedy);
    stats.rootNColsMwis1 += static_cast<int>(pricingSummary.colsPQmwss);
    stats.rootNColsMwis2 += static_cast<int>(pricingSummary.colsPmwss);
    stats.rootNColsExact += static_cast<int>(pricingSummary.colsExact);

    stats.rootTimePool += pricingSummary.timePool;
    stats.rootTimeHeur += pricingSummary.timeGreedy;
    stats.rootTimeMwis1 += pricingSummary.timePQmwss;
    stats.rootTimeMwis2 += pricingSummary.timePmwss;
    stats.rootTimeExact += pricingSummary.timeExact;
  } else {
    stats.otherNodesNCallsPool += static_cast<int>(pricingSummary.callsPool);
    stats.otherNodesNCallsHeur += static_cast<int>(pricingSummary.callsGreedy);
    stats.otherNodesNCallsMwis1 += static_cast<int>(pricingSummary.callsPQmwss);
    stats.otherNodesNCallsMwis2 += static_cast<int>(pricingSummary.callsPmwss);
    stats.otherNodesNCallsExact += static_cast<int>(pricingSummary.callsExact);

    stats.otherNodesNColsPool += static_cast<int>(pricingSummary.colsPool);
    stats.otherNodesNColsHeur += static_cast<int>(pricingSummary.colsGreedy);
    stats.otherNodesNColsMwis1 += static_cast<int>(pricingSummary.colsPQmwss);
    stats.otherNodesNColsMwis2 += static_cast<int>(pricingSummary.colsPmwss);
    stats.otherNodesNColsExact += static_cast<int>(pricingSummary.colsExact);

    stats.otherNodesTimePool += pricingSummary.timePool;
    stats.otherNodesTimeHeur += pricingSummary.timeGreedy;
    stats.otherNodesTimeMwis1 += pricingSummary.timePQmwss;
    stats.otherNodesTimeMwis2 += pricingSummary.timePmwss;
    stats.otherNodesTimeExact += pricingSummary.timeExact;
  }
}

void LP::log_pricing_summary() const {
  if (!params.is_verbose(2)) return;

  size_t totalCalls = pricingSummary.callsPool + pricingSummary.callsGreedy +
                      pricingSummary.callsPQmwss + pricingSummary.callsPmwss +
                      pricingSummary.callsExact;
  size_t totalCols = pricingSummary.colsPool + pricingSummary.colsGreedy +
                     pricingSummary.colsPQmwss + pricingSummary.colsPmwss +
                     pricingSummary.colsExact;
  double totalTime = pricingSummary.timePool + pricingSummary.timeGreedy +
                     pricingSummary.timePQmwss + pricingSummary.timePmwss +
                     pricingSummary.timeExact;

  debugLog << "pricing summary: iters=" << pricingSummary.iters
           << ", total_calls=" << totalCalls << ", total_cols=" << totalCols
           << ", total_time=" << totalTime << std::endl;
  debugLog << "  pool: calls=" << pricingSummary.callsPool
           << ", cols=" << pricingSummary.colsPool
           << ", time=" << pricingSummary.timePool << std::endl;
  debugLog << "  greedy: calls=" << pricingSummary.callsGreedy
           << ", cols=" << pricingSummary.colsGreedy
           << ", time=" << pricingSummary.timeGreedy << std::endl;
  debugLog << "  P,Q-MWSS: calls=" << pricingSummary.callsPQmwss
           << ", cols=" << pricingSummary.colsPQmwss
           << ", time=" << pricingSummary.timePQmwss << std::endl;
  debugLog << "  P-MWSS: calls=" << pricingSummary.callsPmwss
           << ", cols=" << pricingSummary.colsPmwss
           << ", time=" << pricingSummary.timePmwss << std::endl;
  debugLog << "  exact: calls=" << pricingSummary.callsExact
           << ", cols=" << pricingSummary.colsExact
           << ", time=" << pricingSummary.timeExact << std::endl;
}

int LP::pricing_pool(CplexEnv& cenv, IloNumArray& dualsP, IloNumArray& dualsQ) {
  auto startTime = std::chrono::high_resolution_clock::now();
  size_t nPoolCols = 0;
  const size_t k = params.pricingMaxColsPerIter;

  // Min-heap of (cost, pool_index): keeps the k best candidates seen so far.
  // The root is always the worst of the current top-k, so eviction is O(log k).
  using Entry = std::pair<double, size_t>;
  std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> heap;
  for (size_t i = 0; i < pool.size(); ++i) {
    pricingSummary.callsPool++;
    double cost = 0.0;
    for (auto pi : pool[i].ps) cost += dualsP[pi];
    for (auto qj : pool[i].qs) cost -= dualsQ[qj];
    if (cost <= 1 + EPSILON) continue;
    if (heap.size() < k)
      heap.push({cost, i});
    else if (cost > heap.top().first) {
      heap.pop();
      heap.push({cost, i});
    }
  }

  while (!heap.empty()) {
    size_t i = heap.top().second;
    heap.pop();
    add_column(cenv, pool[i]);
    nPoolCols++;
  }

  pricingSummary.colsPool += nPoolCols;
  pricingSummary.timePool += get_elapsed_time(startTime);
  return nPoolCols;
}

int LP::pricing_greedy(CplexEnv& cenv, PricingEnv& penv, IloNumArray& dualsP,
                       IloNumArray& dualsQ, bool enabled) {
  if (!enabled) return 0;
  auto startTime = std::chrono::high_resolution_clock::now();
  size_t nHeurCols = 0;
  const size_t k = params.pricingMaxColsPerIter;

  // Min-heap of Column by cost: keeps the k best candidates seen so far.
  // The root is always the worst of the current top-k, so eviction is O(log k).
  auto cmp = [](const Column& a, const Column& b) { return a.cost > b.cost; };
  std::priority_queue<Column, std::vector<Column>, decltype(cmp)> heap(cmp);
  for (size_t i = 0; i < params.pricingHeur1MaxNCols; ++i) {
    pricingSummary.callsGreedy++;
    auto res = penv.heur_solve(dualsP, dualsQ, params.pricingHeur1Alpha);
    if (res.second != PRICING_STABLE_FOUND) continue;
    if (heap.size() < k)
      heap.push(std::move(res.first));
    else if (res.first.cost > heap.top().cost) {
      heap.pop();
      heap.push(std::move(res.first));
    }
  }

  while (!heap.empty()) {
    Column col = heap.top();
    heap.pop();
    add_column(cenv, col);
    nHeurCols++;
  }

  pricingSummary.colsGreedy += nHeurCols;
  pricingSummary.timeGreedy += get_elapsed_time(startTime);
  return nHeurCols;
}

// P,Q-MWSSP heuristic: the weight of v=(i,j) is \gamma_i - \mu_j
int LP::pricing_P_Q_mwss(CplexEnv& cenv, PricingEnv& penv, IloNumArray& dualsP,
                         IloNumArray& dualsQ, bool enabled) {
  if (!enabled) return 0;
  pricingSummary.callsPQmwss++;
  auto startTime = std::chrono::high_resolution_clock::now();
  size_t nMwis1Cols = 0;
  auto res = penv.mwis_P_Q_solve(dualsP, dualsQ);
  for (auto& [stab, priState] : res) {
    if (priState == PRICING_STABLE_FOUND) {
      add_column(cenv, stab);
      nMwis1Cols++;
    }
  }
  pricingSummary.colsPQmwss += nMwis1Cols;
  pricingSummary.timePQmwss += get_elapsed_time(startTime);
  return nMwis1Cols;
}

// P-MWSSP heuristic: the weight of v=(i,j) is \gamma_i
int LP::pricing_P_mwss(CplexEnv& cenv, PricingEnv& penv, IloNumArray& dualsP,
                       IloNumArray& dualsQ, bool enabled) {
  if (!enabled) return -1;
  pricingSummary.callsPmwss++;
  int ret;
  auto startTime = std::chrono::high_resolution_clock::now();
  auto res = penv.mwis_P_solve(dualsP, dualsQ);
  if (res.second == PRICING_STABLE_FOUND) {
    add_column(cenv, res.first);
    ret = 1;
  } else if (res.second == PRICING_STABLE_NOT_FOUND)
    ret = -1;
  else
    ret = 0;  // Node solved up to optimality

  pricingSummary.colsPmwss += (ret == 1) ? 1 : 0;
  pricingSummary.timePmwss += get_elapsed_time(startTime);
  return ret;
}

int LP::pricing_exact(CplexEnv& cenv, PricingEnv& penv, IloNumArray& dualsP,
                      IloNumArray& dualsQ) {
  pricingSummary.callsExact++;
  auto startTime = std::chrono::high_resolution_clock::now();
  int ret;
  auto res = penv.exact_solve(dualsP, dualsQ);
  // Handle exact outputs
  if (res.second == PRICING_STABLE_FOUND) {
    add_column(cenv, res.first);
    ret = 1;
  } else if (res.second == PRICING_STABLE_NOT_EXIST)
    ret = 0;
  else if (res.second == PRICING_TIME_EXCEEDED)
    ret = -1;
  else if (res.second == PRICING_MEM_EXCEEDED)
    ret = -2;
  else
    ret = -3;

  pricingSummary.colsExact += (ret == 1) ? 1 : 0;
  pricingSummary.timeExact += get_elapsed_time(startTime);
  return ret;
}

void LP::add_column(CplexEnv& cenv, Column& stab) {
  IloNumColumn column = cenv.Xobj(1.0);
  for (auto pi : stab.ps) column += cenv.XrestrP[pi](1.0);
  for (auto qj : stab.qs) {
    if (dpcp.get_Q()[qj].size() > 1) column += cenv.XrestrQ[qj](-1.0);
  }
  cenv.Xvars.add(IloNumVar(column));
  stables.push_back(stab);
  // assert(check_column(stab));
  // print_column(stab);
}

bool LP::check_column(Column& stab) {
  // Check if the column is valid
  for (auto v : stab.stable) {
    if (stab.ps.find(dpcp.get_P_part(v)) == stab.ps.end() ||
        stab.qs.find(dpcp.get_Q_part(v)) == stab.qs.end())
      return false;
  }
  return true;
}

void LP::print_column(Column& stab) {
  std::cout << "Column: [";
  for (auto v : stab.stable) {
    std::cout << " " << dpcp.get_current_id(v) << " = (" << dpcp.get_P_part(v)
              << ", " << dpcp.get_Q_part(v) << ")";
  }
  std::cout << " ]" << std::endl;
  std::cout << "cost: " << stab.cost << std::endl;
  std::cout << "ps: [";
  for (auto pi : stab.ps) std::cout << " " << pi;
  std::cout << " ]" << std::endl;
  std::cout << "qs: [";
  for (auto qj : stab.qs) std::cout << " " << qj;
  std::cout << " ]" << std::endl;
}

// Branching rule of Lucci, Nasini, Tolomei and Torres for PCP
Vertex LP::get_branching_variable_LNTT(const IloNumArray& values) {
  const Vertex null_v = boost::graph_traits<Graph>::null_vertex();
  // Given i in P with index pi,
  // posSnd[pi] = {(v, x): v = (i,j) in V, x = \sum_{S : v in S} x_S }
  std::vector<std::map<Vertex, double>> posSnd(dpcp.get_nP());
  for (auto i : posVars)
    for (auto v : stables[i].stable) {
      size_t pi = dpcp.get_P_part(v);
      if (!posSnd[pi].contains(v)) posSnd[pi][v] = 0.0;
      posSnd[pi][v] += values[i];
    }

  // Branching criterion
  // 1. Let M = min {|posSnd[pi]|: i in P, |posSnd[pi]| > 1}
  // 2. Let P' subset P such that, for all i in P', |posSnd[pi]| = M
  // 3. Choose (v,x) such that
  //    x is the most fractional value among \cup_{i in P'} posSnd[idP_i]
  // 4. Break ties with index of i
  int best_pi = -1;
  Vertex best_v = null_v;
  double best_value = 0.0;
  for (size_t pi = 0; pi < dpcp.get_nP(); ++pi) {
    if (posSnd[pi].size() <= 1)
      continue;
    else if (best_pi == -1 || posSnd[pi].size() < posSnd[best_pi].size()) {
      best_pi = static_cast<int>(pi);
      auto it = find_most_fractional(posSnd[best_pi]);
      best_v = it->first;
      best_value = it->second;
    } else if (posSnd[pi].size() == posSnd[best_pi].size()) {  // Tie
      auto it = find_most_fractional(posSnd[pi]);
      if (std::abs(it->second - 0.5) < std::abs(best_value - 0.5)) {
        best_pi = static_cast<int>(pi);
        best_v = it->first;
        best_value = it->second;
      }
    }
  }

  if (best_v == null_v) {
    // Choose by index
    for (size_t pi = 0; pi < dpcp.get_nP(); ++pi)
      if (dpcp.get_P()[pi].size() > 1) {
        best_v = dpcp.get_P()[pi].front();
        break;
      }
  }

  assert(best_v != null_v);

  // *******************************************************************
  // // Print some statics
  // std::cout << "branching variable: " << best_v << " ["
  //           << dpcp.get_P_part(best_v) << " " <<
  //           dpcp.get_Q_part(best_v)
  //           << "] with size "
  //           << posSnd[dpcp.get_P_part(best_v)].size()
  //           << " and value " << best_value << std::endl
  //           << std::endl;
  // *******************************************************************

  return best_v;
}

// Branching rule of Furuni, Malaguti and Santini for PCP
Vertex LP::get_branching_variable_FMS(const IloNumArray& values) {
  const Vertex null_v = boost::graph_traits<Graph>::null_vertex();
  const Partition& P = dpcp.get_P();
  size_t nP = dpcp.get_nP();

  // Given v in V, stabs[v] is the sum of the values of the variables
  // x_S such that v \in S, i.e., stabs[v] = \sum_{S : v \in S} x_S
  std::map<Vertex, double> stabs;
  for (auto i : posVars)
    for (auto v : stables[i].stable) {
      if (!stabs.contains(v)) stabs[v] = 0.0;
      stabs[v] += values[i];
    }

  // Branching criterion
  // 1. Find pi in P such that P[pi] has the maximum number of partially colored
  //    vertices, i.e, |{v in P[pi]: stabs[v] > 0}| is maximum, as long as this
  //    number is greater than 1. breaking ties by |P[pi]| (the smaller, the
  //    better) and further ties by the index of pi
  size_t best_pi = nP;
  size_t best_nPartial = 0;
  for (size_t pi = 0; pi < nP; ++pi) {
    const auto& Pi = P[pi];
    size_t nPartial = 0;
    for (auto v : Pi) {
      if (stabs.contains(v) && stabs[v] > EPSILON) nPartial++;
    }
    if (nPartial <= 1) continue;
    if (nPartial > best_nPartial) {
      best_nPartial = nPartial;
      best_pi = pi;
    } else if (nPartial == best_nPartial && (Pi.size() < P[best_pi].size())) {
      best_pi = pi;
    }
  }

  // If no such pi exists, choose the first pi such that P[pi] has more than 1
  // vertex
  if (best_pi == nP) {
    // Choose by index
    for (size_t pi = 0; pi < dpcp.get_nP(); ++pi)
      if (dpcp.get_P()[pi].size() > 1) {
        best_pi = pi;
        break;
      }
  }

  // 2. Choose v in P[best_pi] such that stabs[v] is maximum
  Vertex best_v = null_v;
  double best_value = 0.0;
  for (auto v : P[best_pi]) {
    if (stabs.contains(v) && stabs[v] > best_value) {
      best_value = stabs[v];
      best_v = v;
    }
  }

  // If no such v exists, choose the first v in P[best_pi]
  if (best_v == null_v) {
    best_v = P[best_pi].front();
  }

  assert(best_v != null_v);

  // *******************************************************************
  // // Print some statics;
  // std::cout << "branching variable: " << dpcp.get_current_id(best_v) << " ["
  //           << dpcp.get_P_part(best_v) << " " << dpcp.get_Q_part(best_v)
  //           << "] with size " << best_nPartial << " and cardinality "
  //           << dpcp.get_P()[best_pi].size() << " and value " << best_value
  //           << std::endl;
  // *******************************************************************

  return best_v;
}

Col LP::get_lp_solution() {
  assert(state == LP_INTEGER);

  Col originalCol;  // Coloring of the original DPCP instance
  Col currentCol;   // Coloring of the current DPCP instance

  if (dpcp.has_trivial_solution()) {
    Vertex v = *(vertices(dpcp.get_graph()).first);
    size_t v_id = dpcp.get_current_id(v);
    currentCol.set_color(dpcp, v_id, 0);
  } else {
    Color k = 0;
    for (auto i : posVars) {
      for (auto v : stables[i].stable) {
        size_t v_id = dpcp.get_current_id(v);
        // Ignore already colored vertices
        if (currentCol.is_colored(v_id)) continue;
        currentCol.set_color(dpcp, v_id, k);
      }
      ++k;
    }
  }
  assert(currentCol.check_coloring(dpcp));

  // Translate coloring to the original graph
  originalCol = currentCol.translate_coloring(dpcp, origDpcp);

  // Color isolated vertices
  originalCol.color_isolated_vertices(origDpcp, dpcp.get_isolated_vertices());

  assert(originalCol.check_coloring(origDpcp));
  return originalCol;
}

Col LP::get_heur_solution() {
  assert(has_heur_solution());
  Col originalCol;  // Coloring of the original DPCP instance

  // Translate coloring to the original graph
  originalCol = coloring.translate_coloring(dpcp, origDpcp);

  // Color isolated vertices
  originalCol.color_isolated_vertices(origDpcp, dpcp.get_isolated_vertices());

  assert(originalCol.check_coloring(origDpcp));
  return originalCol;
}

void LP::branch(std::vector<std::unique_ptr<LP>>& sons) {
  const Vertex null_v = boost::graph_traits<Graph>::null_vertex();
  Vertex v = branchingVertex;
  assert(v != null_v);

  if (params.is_verbose(2)) {
    debugLog << "Branching on vertex v with id=" << dpcp.get_current_id(v)
             << " (p(v)=" << dpcp.get_P_part(v)
             << ", q(v)=" << dpcp.get_Q_part(v) << ")" << std::endl;
  }

  sons.clear();
  sons.reserve(2);
  sons.emplace_back(std::make_unique<LP>(*this));
  sons.emplace_back(std::make_unique<LP>(*this));

  // *******
  // ** Left branch: v is preselected
  // *******

  DPCPInst& dpcpLeft = sons[0]->get_dpcp_inst();
  Graph& graphLeft = dpcpLeft.get_graph();

  // Map from original vertices to new vertices
  // Necessary to translate the columns of the pool
  std::map<Vertex, Vertex> mapLeft;
  for (auto [v, v_id] : dpcp.get_vertex2CurrentId()) {
    Vertex uLeft = vertex(v_id, graphLeft);
    assert(uLeft != null_v);
    mapLeft.emplace(v, uLeft);
  }

  // Preselect v in the left DPCP instance
  auto itLeftBranch = mapLeft.find(v);
  assert(itLeftBranch != mapLeft.end());
  dpcpLeft.preselect_vertex(itLeftBranch->second);
  if (params.preprocessing) dpcpLeft.preprocess();
  if (params.is_verbose(2)) {
    debugLog << "Left child after preprocessing: |V|="
             << num_vertices(dpcpLeft.get_graph())
             << ", |E|=" << num_edges(dpcpLeft.get_graph())
             << ", |P|=" << dpcpLeft.get_nP() << ", |Q|=" << dpcpLeft.get_nQ()
             << std::endl;
  }

  // *******
  // ** Right branch: v is forbidden
  // *******

  DPCPInst& dpcpRight = sons[1]->get_dpcp_inst();
  Graph& graphRight = dpcpRight.get_graph();

  // Map from original vertices to new vertices
  // Necessary to translate the columns of the pool
  std::map<Vertex, Vertex> mapRight;
  for (auto [v, v_id] : dpcp.get_vertex2CurrentId()) {
    Vertex uRight = vertex(v_id, graphRight);
    assert(uRight != null_v);
    mapRight.emplace(v, uRight);
  }

  // Forbid v in the right DPCP instance
  auto itRightBranch = mapRight.find(v);
  assert(itRightBranch != mapRight.end());
  dpcpRight.forbid_vertex(itRightBranch->second);
  if (params.preprocessing) dpcpRight.preprocess();
  if (params.is_verbose(2)) {
    debugLog << "Right child after preprocessing: |V|="
             << num_vertices(dpcpRight.get_graph())
             << ", |E|=" << num_edges(dpcpRight.get_graph())
             << ", |P|=" << dpcpRight.get_nP() << ", |Q|=" << dpcpRight.get_nQ()
             << std::endl;
  }

  // *******
  // ** Helper lambda to translate columns to a branch
  // *******
  auto translate_columns = [this](const std::vector<size_t>& indices,
                                  const std::map<Vertex, Vertex>& branchMap,
                                  const DPCPInst& branchDpcp) -> Pool {
    Pool result;
    for (size_t i : indices) {
      Column& column = stables[i];
      Column translatedCol;
      for (auto u : column.stable) {
        auto it = branchMap.find(u);
        if (it != branchMap.end()) {
          Vertex uBranch = it->second;
          // uBranch may not be in the branched DPCP instance if it was removed
          // by branching or preprocessing
          if (branchDpcp.has_vertex(uBranch))
            translatedCol.add_vertex(uBranch, branchDpcp.get_P_part(uBranch),
                                     branchDpcp.get_Q_part(uBranch));
        }
      }
      if (translatedCol.stable.size() > 0) result.push_back(translatedCol);
    }
    return result;
  };

  // *******
  // ** Fill pools and late columns
  // *******
  Pool& poolLeft = sons[0]->get_pool();
  Pool& poolRight = sons[1]->get_pool();
  Pool& lateLeft = sons[0]->get_late_columns();
  Pool& lateRight = sons[1]->get_late_columns();

  if (params.inheritColumns > 0) {
    if (params.inheritColumns == 3) {
      // Mode 3: positive columns go directly to LP, non-positive go to pool
      std::vector<size_t> nonPositiveIndices;
      for (size_t i = 0; i < stables.size(); ++i) {
        if (std::find(posVars.begin(), posVars.end(), i) == posVars.end())
          nonPositiveIndices.push_back(i);
      }
      lateLeft = translate_columns(posVars, mapLeft, dpcpLeft);
      poolLeft = translate_columns(nonPositiveIndices, mapLeft, dpcpLeft);
      lateRight = translate_columns(posVars, mapRight, dpcpRight);
      poolRight = translate_columns(nonPositiveIndices, mapRight, dpcpRight);
    } else if (params.inheritColumns == 4) {
      // Mode 4: all columns go directly to LP, pool remains empty
      std::vector<size_t> allIndices;
      allIndices.reserve(stables.size());
      for (size_t i = 0; i < stables.size(); ++i) allIndices.push_back(i);
      lateLeft = translate_columns(allIndices, mapLeft, dpcpLeft);
      lateRight = translate_columns(allIndices, mapRight, dpcpRight);
    } else {
      // Modes 1, 2: columns go to pool (all or positive only)
      std::vector<size_t> inheritIndices;
      if (params.inheritColumns == 1) {
        for (size_t i = 0; i < stables.size(); ++i) inheritIndices.push_back(i);
      } else {
        inheritIndices = posVars;
      }
      poolLeft = translate_columns(inheritIndices, mapLeft, dpcpLeft);
      poolRight = translate_columns(inheritIndices, mapRight, dpcpRight);
    }
  }

  if (params.is_verbose(2)) {
    if (params.inheritColumns == 3 || params.inheritColumns == 4)
      debugLog << "Created child pools: left(pool)=" << poolLeft.size()
               << ", left(late)=" << lateLeft.size()
               << ", right(pool)=" << poolRight.size()
               << ", right(late)=" << lateRight.size() << std::endl;
    else
      debugLog << "Created child pools: left=" << poolLeft.size()
               << ", right=" << poolRight.size() << std::endl;
  }

  return;
}
