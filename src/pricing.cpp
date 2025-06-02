#include "pricing.hpp"
#include "random.hpp"
#include <cfloat>
#include <limits>

// Generic callback: abort when the threshold is exceeded.
// Implementation of early stop.
void ThresholdCallback::check_thresolhd(
    const IloCplex::Callback::Context &context) {
  if (context.isCandidatePoint())
    if (context.getCandidateObjective() > THRESHOLD) {
      // Save stable set
      IloNumArray valY(context.getEnv(), num_vertices(in.graph));
      context.getCandidatePoint(y, valY);
      IloNumArray valWb(context.getEnv(), in.nB);
      context.getCandidatePoint(w, valWb);
      for (auto v : boost::make_iterator_range(vertices(in.graph)))
        if (valY[v] > 0.5) {
          stab.stable.push_back(v);
          stab.as.insert(in.graph[v].first);
        }
      for (size_t iB = 0; iB < in.nB; ++iB)
        if (valWb[iB] > 0.5)
          stab.bs.insert(in.idB2TyB[iB]);
      stab.cost = context.getCandidateObjective();
      // Abort execution
      context.abort();
    }
}

// Implementation of the invoke method.
//
// This is the method that we have to implement to fulfill the
// generic callback contract. CPLEX will call this method during
// the solution process at the places that we asked for.
void ThresholdCallback::invoke(const IloCplex::Callback::Context &context) {
  if (context.inCandidate())
    check_thresolhd(context);
}

PricingEnv::PricingEnv(GraphEnv &in)
    : in(in), stab(), cxenv(), cxmodel(cxenv), y(cxenv, num_vertices(in.graph)),
      w(cxenv, in.nB), cxcons(cxenv), cplex(cxenv), cb(in, stab, y, w),
      contextMask(0), mwis_env(NULL), mwis_pi(NULL), mwis_pi_scalef(1),
      newsets(NULL), nnewsets(0), ecount(0), elist(NULL) {
  exact_init();
  mwis_init();
}

PricingEnv::~PricingEnv() {
  // End CPLEX variables
  cplex.end();
  cxcons.end();
  y.end();
  w.end();
  cxobj.end();
  cxmodel.end();
  cxenv.end();

  // End MWIS variables
  free(elist);
  free(mwis_pi);
  COLORstable_freeenv(&mwis_env);
}

void PricingEnv::exact_init() {

  // Variables
  for (auto v : boost::make_iterator_range(vertices(in.graph))) {
    char name[100];
    snprintf(name, sizeof(name), "y_%ld_%ld", in.graph[v].first,
             in.graph[v].second);
    y[v] = IloBoolVar(cxenv, name);
  }
  for (size_t iB = 0; iB < in.nB; ++iB) {
    char name[100];
    snprintf(name, sizeof(name), "w_%ld", iB);
    w[iB] = IloBoolVar(cxenv, name);
  }

  // Objective
  IloExpr obj(cxenv);
  for (auto v : boost::make_iterator_range(vertices(in.graph)))
    obj += y[v] * 1.0;
  for (size_t iB = 0; iB < in.nB; ++iB)
    obj += w[iB] * (-1.0);
  cxobj = IloMaximize(cxenv, obj);
  cxmodel.add(cxobj);
  obj.end();

  // Constraints
  // (1) \sum_{b \in B: (a,b) \in V} y_a_b <= 1, for all a \in A
  for (size_t ia = 0; ia < in.nA; ++ia) {
    IloExpr restr(cxenv);
    for (size_t v : in.snd[ia])
      restr += y[v];
    cxcons.add(restr <= 1);
    restr.end();
  }

  // (2) y_a_b + y_a'_b' <= 1, for all ((a,b),(a',b')) \in E such that a != a'
  for (auto e : boost::make_iterator_range(edges(in.graph))) {
    auto u = source(e, in.graph);
    auto v = target(e, in.graph);
    if (in.graph[u].first == in.graph[v].first)
      continue;
    IloExpr restr(cxenv);
    restr += y[u] + y[v];
    cxcons.add(restr <= 1);
    restr.end();
  }

  // (3) y_a_b <= w_b, for all (a,b) \in V
  for (auto v : boost::make_iterator_range(vertices(in.graph))) {
    IloExpr restr(cxenv);
    restr += y[v] - w[in.tyB2idB[in.graph[v].second]];
    cxcons.add(restr <= 0);
    restr.end();
  }

  // (4) w_b <= \sum_{v in V^b} y_v, for all b \in B
  for (size_t iB = 0; iB < in.nB; ++iB) {
    IloExpr restr(cxenv);
    restr += w[iB];
    for (auto v : in.fst[iB])
      restr -= y[v];
    cxcons.add(restr <= 0);
    restr.end();
  }

  cxmodel.add(cxcons);

  // Re-export model
  cplex.extract(cxmodel);

  // Set parameters
  cplex.setDefaults();
  cplex.setParam(IloCplex::Param::Parallel, 1); // Deterministic mode
  cplex.setParam(IloCplex::Param::Threads, 1);  // Single thread
  cplex.setOut(cxenv.getNullStream());
  cplex.setWarning(cxenv.getNullStream());

  // Now we get to setting up the generic callback.
  // We instantiate a ThresholdCallback and set the contextMask parameter.
  contextMask |= IloCplex::Callback::Context::Id::Candidate;

  // If contextMask is not zero we add the callback.
  if (contextMask != 0)
    cplex.use(&cb, contextMask);
}

void PricingEnv::mwis_init() {

  // Initalize stable environment
  COLORstable_initenv(&mwis_env, NULL, 0);

  // Intialize vectors of weights
  mwis_pi = (COLORNWT *)COLOR_SAFE_MALLOC(num_vertices(in.graph), COLORNWT);

  // // Initialize maps from original vertices to positive vertices, and
  // viceversa std::iota(toPosVertex.begin(), toPosVertex.end(), 0);
  // std::iota(toVertex.begin(), toVertex.end(), 0);

  // Initialize edge array
  elist = (int *)malloc(sizeof(int) * 2 * num_edges(in.graph));
  for (auto e : boost::make_iterator_range(edges(in.graph))) {
    elist[2 * ecount] = source(e, in.graph);
    elist[2 * ecount++ + 1] = target(e, in.graph);
  }
}

// std::pair<bool, size_t> PricingEnv::get_weights(std::vector<double> &weights,
//                                                 IloNumArray &dualsA,
//                                                 IloNumArray &dualsB) {
//   size_t nPosWeights = 0;
//   bool changed = false;
//   for (auto v : boost::make_iterator_range(vertices(in.graph))) {
//     auto [a, b] = in.graph[v];
//     double weight = dualsA[in.tyA2idA[a]] - dualsB[in.tyB2idB[b]];
//     if ((weight < -PRICING_EPSILON && weights[v] > -PRICING_EPSILON) ||
//         (weight > -PRICING_EPSILON && weights[v] < -PRICING_EPSILON))
//       changed = true;
//     weights[v] = weight;
//     if (weights[v] > -PRICING_EPSILON)
//       nPosWeights++;
//   }
//   return std::make_pair(changed, nPosWeights);
// }

/* Adaptation of  COLOR_double2COLORNWT from exactcolors */
int PricingEnv::double2COLORNWT(COLORNWT nweights[], COLORNWT *scalef,
                                const std::vector<double> &dbl_weights) {
  size_t i;
  double max_dbl_nweight = -DBL_MAX;
  double max_prec_dbl = exp2(DBL_MANT_DIG - 1);
  static const double max_mwiswt = (double)COLORNWT_MAX;
  double dbl_scalef = COLORDBLmin(max_prec_dbl, max_mwiswt);

  dbl_scalef /= (double)dbl_weights.size();

  for (i = 0; i < dbl_weights.size(); ++i) {
    max_dbl_nweight = COLORDBLmax(max_dbl_nweight, dbl_weights[i]);
  }
  dbl_scalef /= COLORDBLmax(1.0, max_dbl_nweight);
  dbl_scalef = floor(dbl_scalef);
  *scalef = (COLORNWT)dbl_scalef;

  for (i = 0; i < dbl_weights.size(); ++i) {
    double weight = dbl_weights[i] * dbl_scalef;
    assert(weight < (double)COLORNWT_MAX);
    nweights[i] = (COLORNWT)weight;
  }
  return 0;
}

// /* Adaptation of  COLOR_double2COLORNWT from exactcolors */
// int LP::double2COLORNWT(COLORNWT nweights[], COLORNWT *scalef,
//                         size_t nPosWeights,
//                         const std::vector<double> &dbl_weights) {
//   size_t i;
//   double max_dbl_nweight = -DBL_MAX;
//   double max_prec_dbl = exp2(DBL_MANT_DIG - 1);
//   static const double max_mwiswt = (double)COLORNWT_MAX;
//   double dbl_scalef = COLORDBLmin(max_prec_dbl, max_mwiswt);

//   // Compute positive dbl weights
//   std::vector<double> dbl_PosWeights;
//   dbl_PosWeights.reserve(nPosWeights);
//   for (auto w : dbl_weights)
//     if (w > -EPSILON)
//       dbl_PosWeights.push_back(w);

//   dbl_scalef /= (double)dbl_PosWeights.size();

//   for (i = 0; i < dbl_PosWeights.size(); ++i) {
//     max_dbl_nweight = COLORDBLmax(max_dbl_nweight, dbl_PosWeights[i]);
//   }
//   dbl_scalef /= COLORDBLmax(1.0, max_dbl_nweight);
//   dbl_scalef = floor(dbl_scalef);
//   *scalef = (COLORNWT)dbl_scalef;

//   for (i = 0; i < dbl_PosWeights.size(); ++i) {
//     double weight = dbl_PosWeights[i] * dbl_scalef;
//     assert(weight < (double)COLORNWT_MAX);
//     nweights[i] = (COLORNWT)weight;
//   }
//   return 0;
// }

std::pair<StableEnv, PRICING_STATE>
PricingEnv::mwis2_solve(IloNumArray &dualsA, IloNumArray &dualsB) {

  // Reset solution counter an cutoff value
  nnewsets = 0;
  mwis_pi_scalef = 10; // high cutoff to force optimality

  // Reset stable
  stab.stable.clear();
  stab.as.clear();
  stab.bs.clear();
  stab.cost = 0.0;

  // Compute vertex weight array
  std::vector<double> weights(num_vertices(in.graph));
  for (auto v : boost::make_iterator_range(vertices(in.graph)))
    weights[v] = dualsA[in.tyA2idA[in.graph[v].first]];

  // Scale weights
  double2COLORNWT(mwis_pi, &mwis_pi_scalef, weights);

  // Solve the MWIS problem
  COLORstable_wrapper(&mwis_env, &newsets, &nnewsets, num_vertices(in.graph),
                      ecount, elist, mwis_pi, mwis_pi_scalef, 0, 0, 2);

  // Recover stable set
  for (int set_i = 0; set_i < nnewsets; ++set_i)
    for (int j = 0; j < newsets[set_i].count; ++j) {
      int v = newsets[set_i].members[j];
      stab.stable.push_back(v);
      auto [a, b] = in.graph[v];
      if (stab.as.insert(a).second)
        stab.cost += dualsA[in.tyA2idA[a]];
      if (stab.bs.insert(b).second)
        stab.cost -= dualsB[in.tyB2idB[b]];
    }

  // Free memory
  free(newsets);

  PRICING_STATE state;
  if (stab.cost > 1 + PRICING_EPSILON)
    state = PRICING_SOLUTION;
  else
    state = PRICING_NO_SOLUTION;

  return std::make_pair(stab, state);
}

// std::pair<StableEnv, PRICING_STATE> PricingEnv::mwis1_solve(IloNumArray
// &dualsA,
//                                                             IloNumArray
//                                                             &dualsB, double
//                                                             timelimit) {
//   // Reset solution counter an cutoff value
//   nnewsets = 0;
//   mwis_pi_scalef = 1;

//   auto [changed, nPosWeights] = get_weights(weights, duals);

//   If the sign of the weight of any vertex has changed...
//   if (changed) {

//     // Reinitialize stable environment
//     COLORstable_freeenv(&mwis_env);
//     COLORstable_initenv(&mwis_env, NULL, 0);

//     // Reindex vertices with positive weight from 0 to nPosWeights-1
//     oldVertices.resize(nPosWeights);
//     size_t i = 0;
//     for (auto v : boost::make_iterator_range(vertices(graph))) {
//       if (weights[v] > -EPSILON) {
//         oldVertices[i] = v;
//         newVertices[v] = i;
//         ++i;
//       } else
//         newVertices[v] = -1;
//     }

//     // Recompute edge list according to the previous indices
//     // Ignore edges that have as endpoint a vertex with a negative weight
//     ecount = 0;
//     for (auto e : boost::make_iterator_range(edges(graph))) {
//       auto v = source(e, graph);
//       auto u = target(e, graph);
//       if (weights[v] < -EPSILON || weights[u] < -EPSILON)
//         continue;
//       elist[2 * ecount] = newVertices[v];
//       elist[2 * ecount++ + 1] = newVertices[u];
//     }
//   }

//   // Scale weights
//   double2COLORNWT(mwis_pi, &mwis_pi_scalef, nPosWeights, weights);
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
// }

std::pair<StableEnv, PRICING_STATE>
PricingEnv::exact_solve(IloNumArray &dualsA, IloNumArray &dualsB) {

  // Update objective coefficients
  IloNumArray y_coefs(cxenv, num_vertices(in.graph));
  for (auto v : boost::make_iterator_range(vertices(in.graph)))
    y_coefs[v] = dualsA[in.tyA2idA[in.graph[v].first]];
  IloNumArray w_coefs(cxenv, in.nB);
  for (size_t iB = 0; iB < in.nB; ++iB)
    w_coefs[iB] = -dualsB[iB];
  cxobj.setLinearCoefs(y, y_coefs);
  cxobj.setLinearCoefs(w, w_coefs);
  cxmodel.add(cxobj);

  // Reset stable
  stab.stable.clear();
  stab.as.clear();
  stab.bs.clear();
  stab.cost = 0.0;

  // Solve
  cplex.extract(cxmodel);
  cplex.setParam(IloCplex::Param::TimeLimit, PRICING_TIMELIMIT);
  cplex.solve();

  // Get final state
  PRICING_STATE state;
  switch (cplex.getCplexStatus()) {
  case IloCplex::CplexStatus::Optimal: {
    // Exit witout abortion means that the optimal value is <= THRESHOLD = 1.1
    // Recover optimal solution
    IloNumArray valY(cxenv, num_vertices(in.graph));
    cplex.getValues(y, valY);
    IloNumArray valW(cxenv, in.nB);
    cplex.getValues(w, valW);
    for (auto v : boost::make_iterator_range(vertices(in.graph)))
      if (valY[v] > 0.5) {
        stab.stable.push_back(v);
        stab.as.insert(in.graph[v].first);
      }
    for (size_t iB = 0; iB < in.nB; ++iB)
      if (valW[iB] > 0.5)
        stab.bs.insert(in.idB2TyB[iB]);
    stab.cost = cplex.getBestObjValue();
    // Classify optimal solution
    if (stab.cost > 1 + PRICING_EPSILON)
      state = PRICING_SOLUTION;
    else
      state = PRICING_NO_SOLUTION;
  } break;
  case IloCplex::CplexStatus::AbortUser:
    // Exit with abortion means that the current value is > THRESHOLD = 1.1
    state = PRICING_SOLUTION;
    // Try to improve the stable set
    for (Vertex v : boost::make_iterator_range(vertices(in.graph))) {
      auto [a, b] = in.graph[v];
      if (dualsA[in.tyA2idA[a]] < PRICING_EPSILON || stab.as.contains(a))
        continue;
      if (!stab.bs.contains(b))
        continue;
      auto it = std::find_if(
          stab.stable.begin(), stab.stable.end(),
          [v, this](auto w) { return edge(v, w, in.graph).second; });
      if (it != stab.stable.end())
        continue;
      stab.stable.push_back(v);
      stab.as.insert(a);
      stab.cost += dualsA[in.tyA2idA[a]];
    }
    break;
  case IloCplex::CplexStatus::AbortTimeLim:
    state = PRICING_TIME_EXCEEDED;
    break;
  case IloCplex::CplexStatus::MemLimFeas:
  case IloCplex::CplexStatus::MemLimInfeas:
    state = PRICING_MEM_EXCEEDED;
    break;
  default:
    state = PRICING_OTHER;
    break;
  }

  return std::make_pair(stab, state);
}

std::pair<StableEnv, PRICING_STATE> PricingEnv::heur_solve(IloNumArray &dualsA,
                                                           IloNumArray &dualsB,
                                                           Vertex start_v) {

  // Reset stable
  stab.stable.clear();
  stab.as.clear();
  stab.bs.clear();
  stab.cost = 0.0;

  // Candidates
  std::list<std::pair<double, Vertex>> candidates;
  TypeA start_a = in.graph[start_v].first;
  TypeA start_b = in.graph[start_v].second;
  for (auto v : boost::make_iterator_range(vertices(in.graph))) {
    auto [a, b] = in.graph[v];
    double cost_a = dualsA[in.tyA2idA[a]];
    double cost_b = dualsB[in.tyB2idB[b]];
    if (v == start_v) {
      stab.stable.push_back(start_v);
      stab.as.insert(start_a);
      stab.bs.insert(start_b);
      stab.cost += cost_a - cost_b;
      continue;
    }
    if (a == start_a || edge(start_v, v, in.graph).second)
      continue;
    if (b == start_b) {
      candidates.push_back(std::make_pair(cost_a, v));
      continue;
    }
    candidates.push_back(std::make_pair(cost_a - cost_b, v));
  }

  while (!candidates.empty()) {

    // Find the best candidate
    auto it_v = std::max_element(candidates.begin(), candidates.end());

    // Discard best candidate with some low probability
    if (rand_double(rng) < 0.05) {
      candidates.erase(it_v);
      continue;
    }

    Vertex v = it_v->second;
    auto [a, b] = in.graph[v];

    // Add the best candidate to the stable set
    stab.stable.push_back(v);
    stab.as.insert(a);
    stab.bs.insert(b);
    stab.cost += it_v->first;

    // Remove and update weights in the candidate list
    for (auto it = candidates.begin(); it != candidates.end();) {
      Vertex u = it->second;
      auto [au, bu] = in.graph[u];
      if (au == a || edge(u, v, in.graph).second) {
        it = candidates.erase(it);
        continue;
      }
      if (bu == b)
        it->first = dualsA[in.tyA2idA[au]];
      ++it;
    }
  }

  PRICING_STATE state;
  if (stab.cost > 1 + PRICING_EPSILON)
    state = PRICING_SOLUTION;
  else
    state = PRICING_NO_SOLUTION;

  return std::make_pair(stab, state);
}