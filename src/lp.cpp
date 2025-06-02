#include "lp.hpp"
#include "pricing.hpp"
#include "random.hpp"

#include <boost/graph/copy.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <cfloat>
#include <limits>
#include <numeric>

LP::LP(const Graph &graph, Pool &_pool, Graph &origGraph, Col *initSol)
    : in(GraphEnv(graph)), stables(), posVars(), objVal(-1.0), initSol(initSol),
      pool(_pool), origGraph(origGraph), vertexMap(), invVertexMap(),
      nTotalPoolCols(0), nTotalHeurCols(0), nTotalMwis2Cols(0),
      nTotalExactCols(0), nRemainingAttemptsPool(MAXFAILSPOOL),
      nRemainingAttemptsHeur(MAXFAILSHEUR),
      nRemainingAttemptsMwis2(MAXFAILSHEUR) {

  // Translate original vertices to current vertices and viceversa
  vertexMap.resize(num_vertices(graph));
  invVertexMap.resize(num_vertices(origGraph), -1);
  for (Vertex v : boost::make_iterator_range(vertices(origGraph))) {
    auto [av, bv] = origGraph[v];
    for (Vertex u : boost::make_iterator_range(vertices(graph))) {
      auto [au, bu] = graph[u];
      if (av == au && bv == bu) {
        vertexMap[u] = v;
        invVertexMap[v] = u;
      }
    }
  }
};

LP::~LP() {
  // for (size_t i = 0; i < stables.size(); ++i) {
  //   if (stables[i]->age >= 0) {
  //     free(stables[i]->members);
  //     free(stables[i]);
  //   } else
  //     delete stables[i];
  // }
}

void LP::add_constraints(CplexEnv &cenv) {
  // Add constraints
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

void LP::add_initial_columns(CplexEnv &cenv) {

  if (initSol != NULL) {
    for (size_t k = 0; k < initSol->get_n_colors(); ++k) {
      StableEnv stab = initSol->get_stable(in.graph, k);
      add_column(cenv, stab);
    }
    return;
  }

  // Add initial columns
  // ¡TODO: this initialization DO NOT work in general cases!
  // Color each b \in B with an unique color
  for (size_t iB = 0; iB < in.nB; ++iB) {
    StableEnv stab;
    stab.stable = in.fst[iB];
    std::set<TypeA> as;
    for (auto v : stab.stable)
      stab.as.insert(in.graph[v].first);
    stab.bs.insert(in.idB2TyB[iB]);
    add_column(cenv, stab);
  }

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
  // *******************************************************************
  // // Print some statics
  // std::cout << "adding column: [";
  // for (auto v : stab.stable) {
  //   auto [a, b] = in.graph[v];
  //   std::cout << " " << v << " = (" << a << ", " << b << ")";
  // }
  // std::cout << " ]" << std::endl;
  // std::cout << "cost: " << stab.cost << std::endl;
  // std::cout << "as: [";
  // for (auto a : stab.as)
  //   std::cout << " " << a;
  // std::cout << " ]" << std::endl;
  // std::cout << "bs: [";
  // for (auto b : stab.bs)
  //   std::cout << " " << b;
  // std::cout << " ]" << std::endl;
  // *******************************************************************
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
  StableEnv newStab;
  for (Vertex v : stab.stable) {
    auto [a, b] = origGraph[v];
    int u = invVertexMap[v];
    if (u < 0)
      return std::make_pair(false, newStab);
    newStab.stable.push_back(u);
    if (newStab.as.insert(a).second)
      newStab.cost += dualsA[in.tyA2idA[a]];
    if (newStab.bs.insert(b).second)
      newStab.cost -= dualsB[in.tyB2idB[b]];
  }
  return std::make_pair(true, newStab);
}

LP_STATE LP::optimize(double timelimit) {

  auto startTime = std::chrono::high_resolution_clock::now();
  LP_STATE state = LP_UNSOLVED;

  // Check if the input is a GCP instance
  if (in.isGCP)
    return solve_GCP(timelimit);

  // Initialize cplex environment
  CplexEnv cenv;
  IloCplex cplex(cenv.Xmodel);
  set_parameters(cenv, cplex);
  add_constraints(cenv);
  add_initial_columns(cenv);

  // Initialize arrays for dual values
  IloNumArray dualsA(cenv.Xenv, in.nA);
  IloNumArray dualsB(cenv.Xenv, in.nB);

  // Initialize pricing environment
  PricingEnv penv(in);

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

    int ret = pricing(cenv, penv, dualsA, dualsB);
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

    values.end();
  }

  cplex.end();
  return state;
}

int LP::pricing(CplexEnv &cenv, PricingEnv &penv, IloNumArray &dualsA,
                IloNumArray &dualsB) {

  // First, look for entering columns in the pool
  if (nRemainingAttemptsPool > 0) {
    size_t nPoolCols = 0;
    for (StableEnv &stab : pool) {
      auto ret = translate_stable_from_pool(stab, dualsA, dualsB);
      if (!ret.first)
        continue;
      if (ret.second.cost > 1 + EPSILON) {
        add_column(cenv, ret.second);
        nPoolCols++;
        nTotalPoolCols++;
        if (nPoolCols > MAXCOLSFROMPOOL)
          break;
      }
    }
    if (nPoolCols > 0) {
      // std::cout << "** Added " << nPoolCols << " columns from pool"
      //           << std::endl;
      return nPoolCols;
    } else
      nRemainingAttemptsPool--;
  }

  std::pair<StableEnv, PRICING_STATE> res;

  // Second, heuristic resolution of pricing
  if (nRemainingAttemptsHeur > 0) {
    size_t nHeurCols = 0;
    for (size_t i = 0; i < MAXCOLSFROMHEUR; ++i) {
      // Random starting vertex for the stable set
      Vertex v = rand_int(rng) % num_vertices(in.graph);
      res = penv.heur_solve(dualsA, dualsB, v);
      if (res.second == PRICING_SOLUTION) {
        add_column(cenv, res.first);
        nHeurCols++;
        nTotalHeurCols++;
      }
    }
    if (nHeurCols > 0) {
      // std::cout << "** Added " << nHeurCols << " columns from heuristic"
      //           << std::endl;
      return nHeurCols;
    } else
      nRemainingAttemptsHeur--;
  }

  // Third, MWSSP heuristic 2: the weight of v is \gamma_v
  if (nRemainingAttemptsMwis2 > 0) {
    res = penv.mwis2_solve(dualsA, dualsB);
    if (res.second == PRICING_SOLUTION) {
      // std::cout << "** Added a column from MWSSP 2" << std::endl;
      add_column(cenv, res.first);
      nTotalMwis2Cols++;
      return 1;
    } else
      nRemainingAttemptsMwis2--;
  }

  // Fourth, exact resolution of pricing
  res = penv.exact_solve(dualsA, dualsB);

  // Handle exact outputs
  if (res.second == PRICING_SOLUTION) {
    // std::cout << "Exact: " << res.first.cost << std::endl;
    add_column(cenv, res.first);
    nTotalExactCols++;
    // std::cout << "** Added a column from ILP" << std::endl;
    return 1;
  } else if (res.second == PRICING_NO_SOLUTION)
    return 0;
  else if (res.second == PRICING_TIME_EXCEEDED)
    return -1;
  else if (res.second == PRICING_MEM_EXCEEDED)
    return -2;
  else
    return -3;
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

size_t LP::get_branching_variable(const IloNumArray &values) {

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
  int best_v = -1;
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

  if (best_v == -1) {
    // Choose by index
    for (size_t iA = 0; iA < in.nA; ++iA)
      if (in.snd[iA].size() > 1) {
        best_v = iA;
        break;
      }
  }

  assert(best_v != -1);

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

void LP::save_solution(Col &col) {
  col.reset_coloring();
  size_t k = 0;
  for (auto i : posVars) {
    for (auto v : stables[i]) {
      // auto [a, b] = in.graph[v];
      // col.set_color(a, b, k);
      col.set_color(vertexMap[v], k);
    }
    ++k;
  }
  assert(col.check_coloring(origGraph));
  return;
}

void LP::branch(std::vector<LP *> &branches) {

  if (N_BRANCHES != 2) {
    // std::cout << "N_BRANCHES != 2: Unimplemented" << std::endl;
    abort();
  }

  // Recover a and b
  auto [a, b] = in.graph[branchVar];

  // *******
  // ** Left branch: (a,b) is colored
  // *******

  // Get the set of vertices to be removed
  // I.e. every vertex (a,b') with b' != b
  std::set<Vertex> removed;
  for (auto v : boost::make_iterator_range(vertices(in.graph)))
    if (in.graph[v].first == a && in.graph[v].second != b)
      removed.insert(v);

  // Create a view of the graph without the vertices to be removed
  using VFilter = boost::is_not_in_subset<std::set<Vertex>>;
  VFilter vFilter(removed);
  using Filter = boost::filtered_graph<Graph, boost::keep_all, VFilter>;
  Filter filtered(in.graph, boost::keep_all(), vFilter);

  // Create a copy of the view (force reindex)
  Graph graph1;
  boost::copy_graph(filtered, graph1);

  // *******
  // ** Right branch: (a,b) is uncolored
  // ** ¡Reuse graph!
  // *******

  Graph graph2(std::move(in.graph));
  clear_vertex(branchVar, graph2);
  remove_vertex(branchVar, graph2);

  // *******
  // ** Create branches
  // *******

  branches.resize(2);
  branches[0] = new LP(graph1, pool, origGraph);
  branches[1] = new LP(graph2, pool, origGraph);

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

/*
void LP::branch(std::vector<LP *> &branches, Vertex v) {

  if (N_BRANCHES != 2) {
    std::cout << "N_BRANCHES != 2: Unimplemented" << std::endl;
    abort();
  }

  // Recover a and b
  auto [a, b] = graph[v];

  // Find some index of interest
  // .. (e-1,u), (e,v1), ..., (e,v), ..., (e,vn), (e+1,w)
  //               i1           i2                  i3
  size_t i1 = 0, i2 = 0, i3 = 0;
  for (; i1 < vcount; ++i1)
    if (cg[vlist[i1]].first == e)
      break;
  for (i2 = i1; i2 < vcount; ++i2)
    if (cg[vlist[i2]].second == v)
      break;
  for (i3 = i2 + 1; i3 < vcount; ++i3)
    if (cg[vlist[i3]].first != e)
      break;

  // *******
  // ** Left branch: e is represented with v
  // *******

  int vcount1 = vcount - i3 + i1 + 1;
  CGVertex *vlist1 = new CGVertex[2 * vcount1];
  memcpy(vlist1, vlist, sizeof(CGVertex) * i1);
  vlist1[i1] = cgv;
  memcpy(vlist1 + i1 + 1, vlist + i3, sizeof(CGVertex) * (vcount - i3));

  int ecount1 = 0;
  int *elist1 = new int[2 * ecount];
  for (int i = 0; i < 2 * ecount; ++i) {
    int u = elist[i];
    if (u >= i1 && u < i3 && u != i2) {
      if (i % 2 == 0)
        ++i;
      continue;
    }
    int u2 = u < i1 ? u : (u == i2 ? i1 : u - i3 + i1 + 1);
    if (i % 2 == 0)
      elist1[2 * ecount1] = u2;
    else {
      elist1[2 * ecount1 + 1] = u2;
      ecount1++;
    }
  }

  LP *lp1 = new LP(hg, cg, vcount1, vlist1, ecount1, elist1, start_t);

  // *******
  // ** Right branch: e is not represented with v
  // *******

  int vcount2 = vcount - 1;
  CGVertex *vlist2 = vlist;
  vlist = NULL;
  memcpy(vlist2 + i2, vlist2 + i2 + 1, sizeof(CGVertex) * (vcount - i2 - 1));

  int ecount2 = 0;
  int *elist2 = new int[2 * ecount];
  for (int i = 0; i < 2 * ecount; ++i) {
    int u = elist[i];
    if (u == i2) {
      if (i % 2 == 0)
        ++i;
      continue;
    }
    int u2 = u < i2 ? u : u - 1;
    if (i % 2 == 0)
      elist2[2 * ecount2] = u2;
    else {
      elist2[2 * ecount2 + 1] = u2;
      ecount2++;
    }
  }

  LP *lp2 = new LP(hg, cg, vcount2, vlist2, ecount2, elist2, start_t);

  branches.resize(2);
  branches[0] = lp1;
  branches[1] = lp2;

  // TODO: Remove
  std::cout << "i1, i2, i3: " << i1 << " " << i2 << " " << i3 << std::endl;

  std::cout << "vcount1: " << vcount1 << std::endl;
  std::cout << "vlist1: ";
  for (int i = 0; i < vcount1; ++i)
    std::cout << i << ":(" << cg[vlist1[i]].first << "," <<
  cg[vlist1[i]].second
              << ") ";
  std::cout << std::endl;

  std::cout << "ecount1: " << ecount1 << std::endl;
  std::cout << "elist1: ";
  for (int i = 0; i < ecount1; ++i)
    std::cout << "(" << elist1[2 * i] << "," << elist1[2 * i + 1] << ") ";
  std::cout << std::endl;

  std::cout << "vcount2: " << vcount2 << std::endl;
  std::cout << "vlist2: ";
  for (int i = 0; i < vcount2; ++i)
    std::cout << i << ":(" << cg[vlist2[i]].first << "," <<
  cg[vlist2[i]].second
              << ") ";
  std::cout << std::endl;

  std::cout << "ecount2: " << ecount2 << std::endl;
  std::cout << "elist2: ";
  for (int i = 0; i < ecount2; ++i)
    std::cout << "(" << elist2[2 * i] << "," << elist2[2 * i + 1] << ") ";
  std::cout << std::endl;

return;
}
*/