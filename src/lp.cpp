#include "lp.hpp"
#include "heur.hpp"
#include "pricing.hpp"
#include "random.hpp"

#include <cfloat>
#include <limits>
#include <numeric>

LP::LP(const Graph &graph, Params &params, Pool &pool, Graph &origGraph,
       bool isRoot)
    : in(GraphEnv(graph, params, isRoot)), params(params), stables(), posVars(),
      objVal(-1.0), initSol(), pool(pool),
      origGraph(origGraph){

          // // Translate original vertices to current vertices and viceversa
          // vertexMap.resize(num_vertices(graph));
          // invVertexMap.resize(num_vertices(origGraph), -1);
          // for (Vertex v : boost::make_iterator_range(vertices(origGraph))) {
          //   auto [av, bv] = origGraph[v];
          //   for (Vertex u : boost::make_iterator_range(vertices(graph))) {
          //     auto [au, bu] = graph[u];
          //     if (av == au && bv == bu) {
          //       vertexMap[u] = v;
          //       invVertexMap[v] = u;
          //     }
          //   }
          // }
      };

LP::~LP() {}

// Add constraints
void LP::add_constraints_and_objective(CplexEnv &cenv) {
  // Add ">= 1" constraints, one for each a \in A
  for (auto i = in.nA; i > 0; --i)
    cenv.XrestrA.add(IloRange(cenv.Xenv, 1.0, IloInfinity));
  // Add "<= 1" constraints, one for each b \in B
  for (auto i = in.nB; i > 0; --i)
    cenv.XrestrB.add(IloRange(cenv.Xenv, -1.0, IloInfinity));
  cenv.Xmodel.add(cenv.XrestrA);
  cenv.Xmodel.add(cenv.XrestrB);
  // Add objective function
  cenv.Xobj = IloMinimize(cenv.Xenv, 0.0);
  cenv.Xmodel.add(cenv.Xobj);
  return;
}

void LP::add_column(CplexEnv &cenv, StableEnv &stab) {
  IloNumColumn column = cenv.Xobj(1.0);
  for (auto a : stab.as)
    column += cenv.XrestrA[in.tyA2idA[a]](1.0);
  for (auto b : stab.bs)
    column += cenv.XrestrB[in.tyB2idB[b]](-1.0);
  cenv.Xvars.add(IloNumVar(column));
  stables.push_back(VertexVector(stab.stable)); // Push a copy
  // print_column(stab);
}

void LP::print_column(StableEnv &stab) {
  std::cout << "Column: [";
  for (auto v : stab.stable) {
    auto [a, b, id] = in.graph[v];
    std::cout << " " << in.getId[v] << " = (" << a << ", " << b << ")";
  }
  std::cout << " ]" << std::endl;
  std::cout << "cost: " << stab.cost << std::endl;
  std::cout << "as: [";
  for (auto a : stab.as)
    std::cout << " " << a;
  std::cout << " ]" << std::endl;
  std::cout << "bs: [";
  for (auto b : stab.bs)
    std::cout << " " << b;
  std::cout << " ]" << std::endl;
}

void LP::add_initial_columns(CplexEnv &cenv) {

  // If there is an initial solution, add its columns
  if (initSol.get_n_colors() > 0) {
    for (size_t k = 0; k < initSol.get_n_colors(); ++k) {
      StableEnv stab = initSol.get_stable(in.graph, k);
      add_column(cenv, stab);
    }
  } else {
    // Add a new dummy variable z with z_a = 1 forall a \in A and z_b = 0 forall
    // b \in B and set a high cost for z
    IloNumColumn z = cenv.Xobj(params.initializationBigWeight);
    for (size_t i = 0; i < in.nA; ++i)
      z += cenv.XrestrA[i](1.0);
    for (size_t i = 0; i < in.nB; ++i)
      z += cenv.XrestrB[i](0.0);
    cenv.Xvars.add(IloNumVar(z));
    stables.push_back(VertexVector()); // Push a null column
  }

  return;
}

void LP::set_parameters(CplexEnv &cenv, IloCplex &cplex) {
  cplex.setDefaults();
  cplex.setParam(IloCplex::Param::Threads, 1);
  cplex.setParam(IloCplex::Param::Parallel, 1);
  cplex.setOut(cenv.Xenv.getNullStream());
  cplex.setWarning(cenv.Xenv.getNullStream());
}

double get_elapsed_time(auto startTime) {
  return std::chrono::duration_cast<std::chrono::seconds>(
             std::chrono::high_resolution_clock::now() - startTime)
      .count();
}

// Translate a stable set from the pool in terms of the vertices of the current
// graph
std::pair<bool, StableEnv> LP::translate_stable_from_pool(StableEnv &stab,
                                                          IloNumArray &dualsA,
                                                          IloNumArray &dualsB) {
  // StableEnv newStab;
  // for (Vertex v : stab.stable) {
  //   auto [a, b] = origGraph[v];
  //   int u = invVertexMap[v];
  //   if (u < 0)
  //     return std::make_pair(false, newStab);
  //   newStab.stable.push_back(u);
  //   if (newStab.as.insert(a).second)
  //     newStab.cost += dualsA[in.tyA2idA[a]];
  //   if (newStab.bs.insert(b).second)
  //     newStab.cost -= dualsB[in.tyB2idB[b]];
  // }
  // return std::make_pair(true, newStab);
  return std::make_pair(false, StableEnv());
}

LP_STATE LP::optimize(double timelimit, Stats &stats) {

  auto startTime = std::chrono::high_resolution_clock::now();
  LP_STATE state = LP_UNSOLVED;

  // Check if the instance is infeasible
  if (in.isInfeasible) {
    std::cout << "Infeasibility detected" << std::endl;
    return LP_INFEASIBLE;
  }

  // Check if the input is a GCP instance
  if (in.isGCP)
    return solve_GCP(timelimit);

  // Initialize cplex environment
  CplexEnv cenv;
  IloCplex cplex(cenv.Xmodel);
  set_parameters(cenv, cplex);
  add_constraints_and_objective(cenv);
  add_initial_columns(cenv);

  // Initialize arrays for dual values
  IloNumArray dualsA(cenv.Xenv, in.nA);
  IloNumArray dualsB(cenv.Xenv, in.nB);

  // Initialize pricing environment
  PricingEnv penv(in, params.pricingExactTimeLimit);

  while (state == LP_UNSOLVED) {

    // Set time limit
    double timelimit2 = timelimit - get_elapsed_time(startTime);
    if (timelimit2 < 0) {
      state = LP_TIME_EXCEEDED;
      break;
    }
    cplex.setParam(IloCplex::Param::TimeLimit, timelimit2);

    // Optimize
    cplex.solve();

    // Handle errores
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
    cplex.getDuals(dualsA, cenv.XrestrA);
    cplex.getDuals(dualsB, cenv.XrestrB);

    // *******************************************************************
    // // Print some statics
    // // Dual values
    // std::cout << "dual values: ";
    // for (size_t i = 0; i < in.nA; ++i)
    //   std::cout << i << "(" << in.idA2TyA[i] << "): " << dualsA[i] << "
    //   ";
    // std::cout << std::endl;
    // for (size_t i = 0; i < in.nB; ++i)
    //   std::cout << i << "(" << in.idB2TyB[i] << "): " << dualsB[i] << "
    //   ";
    // std::cout << std::endl;
    // *******************************************************************

    int ret = pricing(cenv, penv, stats, dualsA, dualsB);
    if (ret > 0)
      continue;
    else if (ret == 0)
      break;
    else if (ret == -1) {
      state = LP_TIME_EXCEEDED_PR;
      break;
    } else if (ret == -2) {
      state = LP_MEM_EXCEEDED_PR;
      break;
    }
  }

  if (state == LP_UNSOLVED) {

    // Recover primal values and objective value
    IloNumArray values = IloNumArray(cenv.Xenv, cenv.Xvars.getSize());
    cplex.getValues(values, cenv.Xvars);
    objVal = cplex.getObjValue();

    // *******************************************************************
    // // Print some statics
    // // Primal values
    // std::cout << "primal values: ";
    // for (int i = 0; i < cenv.Xvars.getSize(); ++i)
    //   std::cout << values[i] << " ";
    // std::cout << std::endl;
    // // Objective function
    // std::cout << "objective value: " << objVal << std::endl;
    // // Number of columns added
    // std::cout << "Columnas del pool: " << nPoolColsTotal << std::endl;
    // std::cout << "Columnas heurísticas: " << nHeurColsTotal << std::endl;
    // std::cout << "Columnas exactas: " << nExactColsTotal << std::endl;
    // *******************************************************************

    // Check if dummy column is not in the basis
    // We assume that the dummy column is the first one added
    if (values[0] < EPSILON) {

      // Integrality check
      state = LP_INTEGER;
      for (int i = 0; i < cenv.Xvars.getSize(); ++i) {
        if (values[i] < EPSILON)
          continue;
        else if (values[i] < 1 - EPSILON)
          state = LP_FRACTIONAL;
        posVars.push_back(i);
      }

      if (state == LP_FRACTIONAL) {
        // Find branching variable
        branchVar = get_branching_variable(values);
      } else if (state == LP_INTEGER)
        objVal = posVars.size();
    }
    // If the dummy columns is in the basis and the objective value is greater
    // than W, then the instance is infeasible
    else if (objVal > std::min(in.nA, in.nB) + EPSILON)
      state = LP_INFEASIBLE;
    // Otherwise, the initialization failed
    else {
      std::cout << "Initialization failed" << std::endl;
      abort();
    }
    values.end();
  }

  cplex.end();
  return state;
}

int LP::pricing_pool(CplexEnv &cenv, Stats &stats, IloNumArray &dualsA,
                     IloNumArray &dualsB) {
  if (!params.usePool)
    return 0;
  size_t nPoolCols = 0;
  for (StableEnv &stab : pool) {
    auto ret = translate_stable_from_pool(stab, dualsA, dualsB);
    // Can the stable set be used in the current graph?
    if (!ret.first)
      continue;
    // Is the stable set entering?
    if (ret.second.cost > 1 + EPSILON) {
      add_column(cenv, ret.second);
      nPoolCols++;
      stats.nColsPool++;
      // if (nPoolCols > MAXCOLS)
      //   break;
    }
  }
  return nPoolCols;
}

int LP::pricing_greedy(CplexEnv &cenv, PricingEnv &penv, Stats &stats,
                       IloNumArray &dualsA, IloNumArray &dualsB) {
  if (!params.pricingHeur1)
    return 0;
  size_t nHeurCols = 0;
  for (size_t i = 0; i < params.pricingHeur1MaxNCols; ++i) {
    auto res = penv.heur_solve(dualsA, dualsB);
    stats.nCallsHeur++;
    if (res.second == PRICING_STABLE_FOUND) {
      add_column(cenv, res.first);
      nHeurCols++;
      stats.nColsHeur++;
    }
  }
  return nHeurCols;
}

// MWSSP heuristic I: the weight of (a,b) is \gamma_a - \mu_b
int LP::pricing_mwss1(CplexEnv &cenv, PricingEnv &penv, Stats &stats,
                      IloNumArray &dualsA, IloNumArray &dualsB) {
  if (!params.pricingHeur2)
    return 0;
  size_t nMwis1Cols = 0;
  auto res = penv.mwis1_solve(dualsA, dualsB);
  stats.nCallsMWis1++;
  for (auto &[stab, state] : res) {
    if (state == PRICING_STABLE_FOUND) {
      add_column(cenv, stab);
      nMwis1Cols++;
      stats.nColsMwis1++;
    }
  }
  return nMwis1Cols;
}

// MWSSP heuristic II: the weight of (a,b) is \gamma_a
int LP::pricing_mwss2(CplexEnv &cenv, PricingEnv &penv, Stats &stats,
                      IloNumArray &dualsA, IloNumArray &dualsB) {
  if (!params.pricingHeur3)
    return 0;
  auto res = penv.mwis2_solve(dualsA, dualsB);
  stats.nCallsMWis2++;
  if (res.second == PRICING_STABLE_FOUND) {
    add_column(cenv, res.first);
    stats.nColsMwis2++;
    return 1;
  } else if (res.second == PRICING_STABLE_NOT_FOUND)
    return -1;
  else
    return 0; // Node solved up to optimality
}

int LP::pricing_exact(CplexEnv &cenv, PricingEnv &penv, Stats &stats,
                      IloNumArray &dualsA, IloNumArray &dualsB) {
  auto res = penv.exact_solve(dualsA, dualsB);
  stats.nCallsExact++;

  // Handle exact outputs
  if (res.second == PRICING_STABLE_FOUND) {
    add_column(cenv, res.first);
    stats.nColsExact++;
    return 1;
  } else if (res.second == PRICING_STABLE_NOT_EXIST)
    return 0;
  else if (res.second == PRICING_TIME_EXCEEDED)
    return -1;
  else if (res.second == PRICING_MEM_EXCEEDED)
    return -2;
  else
    return -3;
}

int LP::pricing(CplexEnv &cenv, PricingEnv &penv, Stats &stats,
                IloNumArray &dualsA, IloNumArray &dualsB) {

  int addedCols = 0;

  // First, look for entering columns in the pool
  addedCols = pricing_pool(cenv, stats, dualsA, dualsB);
  if (addedCols > 0)
    return addedCols;

  // Second, look for entering columns with a greedy heuristic
  addedCols = pricing_greedy(cenv, penv, stats, dualsA, dualsB);
  if (addedCols > 0)
    return addedCols;

  // Third, MWSSP heuristic I
  addedCols = pricing_mwss1(cenv, penv, stats, dualsA, dualsB);
  if (addedCols > 0)
    return addedCols;

  // Fourth, MWSSP heuristic II
  addedCols = pricing_mwss2(cenv, penv, stats, dualsA, dualsB);
  if (addedCols >= 0)
    return addedCols;

  // Fifth, exact resolution of pricing
  return pricing_exact(cenv, penv, stats, dualsA, dualsB);
}

/* Solve a graph coloring problem instance with exactcolors */
LP_STATE LP::solve_GCP(double timelimit) {
  LP_STATE state;

  COLORproblem colorproblem;
  COLORparms *parms = &(colorproblem.parms);
  colordata *cd = &(colorproblem.root_cd);
  int ncolors = 0;

  if (timelimit < 0)
    return LP_TIME_EXCEEDED;

  // Build GCP graph
  GCPGraph graphGCP;
  get_gcp_graph(in.graph, graphGCP, in.tyB2idB, in.idB2TyB);

  // Get edge list
  int ecount = 0;
  int elist[2 * num_edges(graphGCP)];
  for (auto e : boost::make_iterator_range(edges(graphGCP))) {
    elist[2 * ecount] = source(e, graphGCP);
    elist[2 * ecount++ + 1] = target(e, graphGCP);
  }

  // Initialization
  COLORset *colorclasses = (COLORset *)NULL;
  COLORproblem_init_with_graph(&colorproblem, num_vertices(graphGCP), ecount,
                               elist);
  cd->id = 0;
  colorproblem.ncolordata = 1;
  parms->branching_cpu_limit = timelimit;

  int COLORproblem_init_with_graph(COLORproblem * problem, int ncount,
                                   int ecount, const int elist[]);

  // Find exact coloring
  COLORexact_coloring(&colorproblem, &ncolors, &colorclasses);

  // Optimality check
  if (cd->lower_bound == cd->upper_bound) {
    state = LP_INTEGER;
    objVal = ncolors;
    // Save solution
    // Stable sets need to be translated into vertices of the original graph
    for (int i = 0; i < ncolors; ++i) {
      colorclasses[i].age = 0;

      // Find the correct size
      int newCount = 0;
      for (int j = 0; j < colorclasses[i].count; ++j)
        newCount += in.fst[colorclasses[i].members[j]].size();

      VertexVector stable;
      stable.reserve(newCount);

      // // Alloc memory for the translated stable set
      // COLORset *nset = (COLORset *)malloc(sizeof(COLORset));
      // memcpy(nset, &colorclasses[i], sizeof(COLORset));
      // nset->members = (int *)malloc(sizeof(int) * newCount);
      // nset->count = newCount;

      // Write translated stable set
      for (int j = 0; j < colorclasses[i].count; ++j)
        for (auto v : in.fst[colorclasses[i].members[j]])
          // nset->members[l++] = v;
          stable.push_back(v);

      // Add the translated stable set
      stables.push_back(stable);
      // stables.push_back(nset);
      posVars.push_back(i);
    }
  } else
    state = LP_TIME_EXCEEDED;

  COLORproblem_free(&colorproblem);
  COLORfree_sets(&colorclasses, &ncolors);
  COLORlp_free_env();
  return state;
}

auto find_most_fractional(std::map<Vertex, double> &m) {
  return std::max_element(m.begin(), m.end(),
                          [](const std::pair<Vertex, double> &a,
                             const std::pair<Vertex, double> &b) -> bool {
                            return std::abs(a.second - 0.5) >
                                   std::abs(b.second - 0.5);
                          });
}

Vertex LP::get_branching_variable(const IloNumArray &values) {

  // Given a in A with index i_a,
  // posSnd[i_a] = {(v, x): v = (a,b) \in V, x = \sum_{S : v \in S} x_S }
  std::vector<std::map<Vertex, double>> posSnd(in.nA);
  for (auto i : posVars)
    for (auto v : stables[i]) {
      size_t iA = in.tyA2idA[in.graph[v].first];
      if (!posSnd[iA].contains(v))
        posSnd[iA][v] = 0.0;
      posSnd[iA][v] += values[i];
    }

  // Branching criterion
  // 1. Let M = min {|posSnd[i_a]|: a \in A, |posSnd[i_a]| > 1}
  // 2. Let A' \subset A such that, for all a \in A', |posSnd[i_a]| = M
  // 3. Choose (v,x) such that
  //    x is the most fracional value among \cup_{a \in A'} posSnd[i_a]
  // 4. Break ties with index of a
  int best_a = -1;
  Vertex best_v = NULL;
  double best_value = 0.0;
  for (size_t i = 0; i < in.nA; ++i) {
    if (posSnd[i].size() <= 1)
      continue;
    else if (best_a == -1 || posSnd[i].size() < posSnd[best_a].size()) {
      best_a = i;
      auto it = find_most_fractional(posSnd[best_a]);
      best_v = it->first;
      best_value = it->second;
    } else if (posSnd[i].size() == posSnd[best_a].size()) { // Tie
      auto it = find_most_fractional(posSnd[i]);
      if (std::abs(it->second - 0.5) < std::abs(best_value - 0.5)) {
        best_a = i;
        best_v = it->first;
        best_value = it->second;
      }
    }
  }

  if (best_v == NULL) {
    // Choose by index
    for (size_t iA = 0; iA < in.nA; ++iA)
      if (in.snd[iA].size() > 1) {
        best_v = in.snd[iA].front();
        break;
      }
  }

  assert(best_v != NULL);

  // *******************************************************************
  // // Print some statics
  // std::cout << "branching variable: " << best_v << " ["
  //           << in.graph[best_v].first << " " << in.graph[best_v].second
  //           << "] with size "
  //           << posSnd[in.tyA2idA[in.graph[best_v].first]].size()
  //           << " and value " << best_value << std::endl
  //           << std::endl;
  // *******************************************************************

  return best_v;
}

bool LP::solve_heur(double &obj_value, double &time, bool isRoot) {
  if (!params.initializationUseHeur)
    return false;
  Stats stats;
  if (isRoot)
    stats = dpcp_2_step_semigreedy_heur(in, initSol, 100);
  else
    stats = dpcp_2_step_greedy_heur(in, initSol);
  time = stats.time;
  if (stats.state == FEASIBLE) {
    obj_value = initSol.get_n_colors();
    return true;
  }
  return false;
}

void LP::save_lp_solution(Col &col) {

  // Reset coloring
  col.reset_coloring();

  // Recover the coloring of the current graph
  Col currentCol;
  Color k = 0;
  for (auto i : posVars) {
    for (auto v : stables[i])
      currentCol.set_color(in.graph, v, k);
    ++k;
  }
  assert(currentCol.check_coloring(in.graph));

  // Translate coloring to the original graph
  currentCol.translate_coloring(in.graph, origGraph, col);

  // Color isolated vertices
  currentCol.color_isolated_vertices(in.isolated, col, origGraph);

  assert(col.check_coloring(origGraph));
  return;
}

void LP::save_heur_solution(Col &col) {
  col.reset_coloring();
  col = initSol;
  return;
}

void LP::branch(std::vector<LP *> &branches) {

  // *******
  // ** Left branch: v is colored
  // *******

  // Create a copy of the graph
  Graph gcopy = graph_copy(in.graph, in.getId);
  Vertex branchVarCopy = vertex(in.getId[branchVar], gcopy);
  vertex_branching1(gcopy, branchVarCopy);

  // *******
  // ** Right branch: v is uncolored
  // ** ¡Reuse graph!
  // *******

  vertex_branching2(in.graph, branchVar);

  // *******
  // ** Create branches
  // *******

  branches.resize(2);
  branches[0] = new LP(gcopy, params, pool, origGraph);
  branches[1] = new LP(in.graph, params, pool, origGraph);

  // for (auto v : boost::make_iterator_range(vertices(graph1)))
  //   std::cout << v << ": (" << graph1[v].first << "," << graph1[v].second
  //             << ") ";
  // std::cout << std::endl;

  // for (auto v : boost::make_iterator_range(vertices(graph2)))
  //   std::cout << v << ": (" << graph2[v].first << "," << graph2[v].second
  //             << ") ";
  // std::cout << std::endl;

  return;
}