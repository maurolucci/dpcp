#include "lp.hpp"
#include "pricing.hpp"

#include <boost/graph/copy.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <cfloat>
#include <limits>
#include <numeric>

LP::LP(const Graph &graph) : LP(Graph{graph}){};

LP::LP(const Graph &&graph)
    : in(GraphEnv(graph)), stables(), posVars(), objVal(-1.0){};

LP::~LP() {
  // for (size_t i = 0; i < stables.size(); ++i) {
  //   if (stables[i]->age >= 0) {
  //     free(stables[i]->members);
  //     free(stables[i]);
  //   } else
  //     delete stables[i];
  // }
}

void LP::initialize(CplexEnv &cenv) {

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

  // Add initial columns
  // ¡TODO: this initialization DO NOT work in general cases!
  // Color each b \in B with an unique color
  for (size_t iB = 0; iB < in.nB; ++iB) {
    StableEnv stab;
    stab.stable = in.fst[iB];
    std::set<TypeA> as;
    for (auto v : stab.stable)
      as.insert(in.graph[v].first);
    for (auto a : as)
      stab.as.push_back(a);
    stab.bs.push_back(in.idB2TyB[iB]);
    add_column(cenv, stab);
  }

  return;
}

// void LP::add_column(CplexEnv &cenv, COLORset *newset) {
//   std::vector<bool> cA(in.nA);
//   std::vector<bool> cB(in.nB);
//   for (int i = 0; i < newset->count; ++i) {
//     auto [a, b] = in.graph[newset->members[i]];
//     cA[in.tyA2idA[a]] = true;
//     cB[in.tyB2idB[b]] = true;
//   }
//   IloNumColumn column = cenv.Xobj(1.0);
//   for (size_t i = 0; i < cA.size(); ++i)
//     if (cA[i])
//       column += cenv.Xrestr[i](1.0);
//   for (size_t i = 0; i < cB.size(); ++i)
//     if (cB[i])
//       column += cenv.Xrestr[cA.size() + i](-1.0);
//   cenv.Xvars.add(IloNumVar(column));
//   stables.push_back(newset);

//   // *******************************************************************
//   // Print some statics
//   std::cout << "adding column: [";
//   for (size_t i = 0; i < cA.size(); ++i)
//     if (cA[i])
//       std::cout << " " << in.idA2TyA[i];
//   std::cout << " ] [";
//   for (size_t i = 0; i < cB.size(); ++i)
//     if (cB[i])
//       std::cout << " " << in.idB2TyB[i];
//   std::cout << " ]" << std::endl;
//   // *******************************************************************
//   std::cout << "adding column: [";
//   for (int i = 0; i < newset->count; ++i) {
//     auto [a, b] = in.graph[newset->members[i]];
//     std::cout << " " << newset->members[i] << " = (" << a << ", " << b <<
//     ")";
//   }
//   std::cout << " ]" << std::endl;
//   // *******************************************************************
// }

void LP::add_column(CplexEnv &cenv, StableEnv &stab) {
  IloNumColumn column = cenv.Xobj(1.0);
  for (auto a : stab.as)
    column += cenv.XrestrA[in.tyA2idA[a]](1.0);
  for (auto b : stab.bs)
    column += cenv.XrestrB[in.tyB2idB[b]](-1.0);
  cenv.Xvars.add(IloNumVar(column));
  stables.push_back(VertexVector(stab.stable)); // Push a copy
  // *******************************************************************
  // Print some statics
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

std::pair<bool, size_t> LP::get_weights(std::vector<double> &weights,
                                        IloNumArray &duals) {
  size_t nPosWeights = 0;
  bool changed = false;
  for (auto v : boost::make_iterator_range(vertices(in.graph))) {
    auto [a, b] = in.graph[v];
    double weight = duals[in.tyA2idA[a]] - duals[in.nA + in.tyB2idB[b]];
    if ((weight < -EPSILON && weights[v] > -EPSILON) ||
        (weight > -EPSILON && weights[v] < -EPSILON))
      changed = true;
    weights[v] = weight;
    if (weights[v] > -EPSILON)
      nPosWeights++;
  }
  return std::make_pair(changed, nPosWeights);
}

LP_STATE LP::optimize(double timelimit) {

  auto startTime = std::chrono::high_resolution_clock::now();
  LP_STATE state = LP_UNSOLVED;

  // MWISenv *mwis_env = NULL;
  // COLORNWT *mwis_pi = NULL;
  // COLORNWT mwis_pi_scalef = 1;

  // COLORset *newsets = NULL;
  // int nnewsets = 0;

  // Check if the input is a GCP instance
  if (in.isGCP)
    return solve_GCP(timelimit);

  // Initialize cplex environment
  CplexEnv cenv;
  IloCplex cplex(cenv.Xmodel);
  set_parameters(cenv, cplex);
  initialize(cenv);
  IloNumArray dualsA(cenv.Xenv, in.nA);
  IloNumArray dualsB(cenv.Xenv, in.nB);

  // Initialize pricing environment
  PricingEnv penv(in);

  // // Initalize stable environment
  // COLORstable_initenv(&mwis_env, NULL, 0);

  // // Intialize vectors of weights
  // mwis_pi = (COLORNWT *)COLOR_SAFE_MALLOC(num_vertices(graph), COLORNWT);
  // std::vector<double> weights(num_vertices(graph), 0.0);

  // // Map from new vertices (with positive weight) to original vertices,
  // // and viceversa
  // std::vector<int> newVertices(num_vertices(graph));
  // std::iota(newVertices.begin(), newVertices.end(), 0);
  // std::vector<Vertex> oldVertices(num_vertices(graph));
  // std::iota(oldVertices.begin(), oldVertices.end(), 0);

  // // Initialize edge array
  // int ecount = 0;
  // int elist[2 * num_edges(graph)];
  // for (auto e : boost::make_iterator_range(edges(graph))) {
  //   elist[2 * ecount] = source(e, graph);
  //   elist[2 * ecount++ + 1] = target(e, graph);
  // }

  while (state == LP_UNSOLVED) {

    // // Reset solution counter an cutoff value
    // nnewsets = 0;
    // mwis_pi_scalef = 1;

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
    for (size_t iB = 0; iB < in.nB; ++iB)
      dualsB[iB] *= -1.0;

    // // Compute vertex weight array
    // auto [changed, nPosWeights] = get_weights(weights, duals);

    // // If the sign of the weight of any vertex has changed...
    // if (changed) {

    //   // Reinitialize stable environment
    //   COLORstable_freeenv(&mwis_env);
    //   COLORstable_initenv(&mwis_env, NULL, 0);

    //   // Reindex vertices with positive weight from 0 to nPosWeights-1
    //   oldVertices.resize(nPosWeights);
    //   size_t i = 0;
    //   for (auto v : boost::make_iterator_range(vertices(graph))) {
    //     if (weights[v] > -EPSILON) {
    //       oldVertices[i] = v;
    //       newVertices[v] = i;
    //       ++i;
    //     } else
    //       newVertices[v] = -1;
    //   }

    //   // Recompute edge list according to the previous indices
    //   // Ignore edges that have as endpoint a vertex with a negative weight
    //   ecount = 0;
    //   for (auto e : boost::make_iterator_range(edges(graph))) {
    //     auto v = source(e, graph);
    //     auto u = target(e, graph);
    //     if (weights[v] < -EPSILON || weights[u] < -EPSILON)
    //       continue;
    //     elist[2 * ecount] = newVertices[v];
    //     elist[2 * ecount++ + 1] = newVertices[u];
    //   }
    // }

    // // Scale weights
    // double2COLORNWT(mwis_pi, &mwis_pi_scalef, nPosWeights, weights);

    // *******************************************************************
    // // Print some statics
    // // Dual values
    // std::cout << "dual values: ";
    // for (size_t i = 0; i < in.nA; ++i)
    //   std::cout << i << "(" << in.idA2TyA[i] << "): " << dualsA[i] << " ";
    // std::cout << std::endl;
    // for (size_t i = 0; i < in.nB; ++i)
    //   std::cout << i << "(" << in.idB2TyB[i] << "): " << dualsB[i] << " ";
    // std::cout << std::endl;
    // // Weights
    // std::cout << "weights: ";
    // for (size_t i = 0; i < weights.size(); ++i)
    //   std::cout << weights[i] << " ";
    // std::cout << std::endl;
    // // Positive weights
    // std::cout << "number of positive weights: " << nPosWeights << " of "
    //           << num_vertices(graph) << std::endl;
    // // Scaled weights
    // std::cout << "mwis_pi: ";
    // for (size_t i = 0; i < nPosWeights; ++i)
    //   std::cout << mwis_pi[i] << " ";
    // std::cout << std::endl;
    // std::cout << "mwis_pi_scalef: " << mwis_pi_scalef << std::endl;
    // // Edge array
    // std::cout << "edges: ";
    // for (int i = 0; i < ecount; ++i)
    //   std::cout << "(" << elist[2 * i] << "," << elist[2 * i + 1] << ") ";
    // std::cout << std::endl;
    // // Changed
    // std::cout << "changed: " << changed << std::endl;
    // *******************************************************************

    // if (ecount == 0) {
    //   // Edge-less graphs raise error in COLORstable_wrapper
    //   // So, they are manually solved
    //   nnewsets = 1;
    //   newsets = (COLORset *)malloc(sizeof(COLORset) * nnewsets);
    //   newsets[0].next = NULL;
    //   newsets[0].count = nPosWeights;
    //   newsets[0].members = (int *)malloc(sizeof(int) * newsets[0].count);
    //   for (int i = 0; i < newsets[0].count; ++i)
    //     newsets[0].members[i] = i;
    // } else
    //   // Solve the MWIS problem
    //   COLORstable_wrapper(&mwis_env, &newsets, &nnewsets, nPosWeights,
    //   ecount,
    //                       elist, mwis_pi, mwis_pi_scalef, 0, 0, 2);

    // // Add column/s
    // for (int set_i = 0; set_i < nnewsets; ++set_i) {
    //   // Translate stable set into old vertices
    //   for (int j = 0; j < newsets[set_i].count; ++j) {
    //     newsets[set_i].members[j] = oldVertices[newsets[set_i].members[j]];
    //   }
    //   newsets[set_i].age = 0;
    //   // Local copy of the translated stable set
    //   COLORset *nset = (COLORset *)malloc(sizeof(COLORset));
    //   memcpy(nset, &newsets[set_i], sizeof(COLORset));
    //   // Add column
    //   add_column(cenv, nset);
    // }

    // // Free memory
    // free(newsets);

    std::pair<StableEnv, PRICING_STATE> res;

    // First, heuristic resolution of pricing
    penv.heur_init(dualsA, dualsB);
    size_t added = 0;
    for (size_t iA = 0; iA < in.nA; ++iA) {
      res = penv.heur_solve(dualsA, dualsB, in.idA2TyA[iA]);
      // Non-optimality check
      if (res.first.cost > 1 + EPSILON) {
        add_column(cenv, res.first);
        added++;
      }
    }

    if (added > 0) {
      // std::cout << "Heuristic added " << added << " columns" << std::endl;
      continue;
    }

    // Second, exact resolution of pricing
    res = penv.exact_solve(dualsA, dualsB, timelimit2);

    // Handle errors
    if (res.second == PRICING_TIME_EXCEEDED) {
      state = LP_TIME_EXCEEDED;
      break;
    } else if (res.second == PRICING_MEM_EXCEEDED) {
      state = LP_MEM_EXCEEDED;
      break;
    }

    // Non-optimality check
    if (res.first.cost > 1 + EPSILON) {
      // std::cout << "Exact" << std::endl;
      add_column(cenv, res.first);
      continue;
    }

    // Optimality proved
    break;
  }

  if (state != LP_TIME_EXCEEDED && state != LP_MEM_EXCEEDED) {

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

  // free(mwis_pi);
  // COLORstable_freeenv(&mwis_env);
  cplex.end();
  return state;
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

/* Adaptation of  COLOR_double2COLORNWT from exactcolors */
int LP::double2COLORNWT(COLORNWT nweights[], COLORNWT *scalef,
                        size_t nPosWeights,
                        const std::vector<double> &dbl_weights) {
  size_t i;
  double max_dbl_nweight = -DBL_MAX;
  double max_prec_dbl = exp2(DBL_MANT_DIG - 1);
  static const double max_mwiswt = (double)COLORNWT_MAX;
  double dbl_scalef = COLORDBLmin(max_prec_dbl, max_mwiswt);

  // Compute positive dbl weights
  std::vector<double> dbl_PosWeights;
  dbl_PosWeights.reserve(nPosWeights);
  for (auto w : dbl_weights)
    if (w > -EPSILON)
      dbl_PosWeights.push_back(w);

  dbl_scalef /= (double)dbl_PosWeights.size();

  for (i = 0; i < dbl_PosWeights.size(); ++i) {
    max_dbl_nweight = COLORDBLmax(max_dbl_nweight, dbl_PosWeights[i]);
  }
  dbl_scalef /= COLORDBLmax(1.0, max_dbl_nweight);
  dbl_scalef = floor(dbl_scalef);
  *scalef = (COLORNWT)dbl_scalef;

  for (i = 0; i < dbl_PosWeights.size(); ++i) {
    double weight = dbl_PosWeights[i] * dbl_scalef;
    assert(weight < (double)COLORNWT_MAX);
    nweights[i] = (COLORNWT)weight;
  }
  return 0;
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
  // std::cout << "branching variable: " << best_v << " [" <<
  // graph[best_v].first
  //           << " " << graph[best_v].second << "] with size "
  //           << posSnd[tyA2idA[graph[best_v].first]].size() << " and value "
  //           << best_value << std::endl
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
      col.set_color(v, k);
    }
    ++k;
  }
  assert(col.check_coloring(in.graph));
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
  branches[0] = new LP(graph1);
  branches[1] = new LP(graph2);

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