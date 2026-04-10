#include "pricing.hpp"

#include <cfloat>
#include <limits>

#include "random.hpp"

// Generic callback: abort when the threshold is exceeded.
// Implementation of early stop.
void ThresholdCallback::check_threshold(
    const IloCplex::Callback::Context& context) {
  if (context.isCandidatePoint())
    if (context.getCandidateObjective() > THRESHOLD) {
      // Save stable set
      IloNumArray valY(context.getEnv(), num_vertices(dpcp.get_graph()));
      context.getCandidatePoint(y, valY);
      IloNumArray valWb(context.getEnv(), dpcp.get_nQ());
      context.getCandidatePoint(w, valWb);
      for (size_t qj = 0; qj < dpcp.get_nQ(); ++qj) {
        if (valWb[qj] < 0.5) continue;
        // Postprocessing: force w_j <= \sum_{v in Q[j]} y_v
        size_t sum_y = 0;
        for (auto v : dpcp.get_Q()[qj])
          if (valY[dpcp.get_current_id(v)] > 0.5) sum_y += 1;
        if (sum_y == 0)
          valWb[qj] = 0.0;
        else
          stab.qs.insert(qj);
      }
      for (auto v : boost::make_iterator_range(vertices(dpcp.get_graph())))
        if (valY[dpcp.get_current_id(v)] > 0.5) {
          stab.stable.push_back(v);
          stab.ps.insert(dpcp.get_P_part(v));
        }
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
void ThresholdCallback::invoke(const IloCplex::Callback::Context& context) {
  if (context.inCandidate()) check_threshold(context);
}

PricingEnv::PricingEnv(DPCPInst& dpcpRef, double exactTimeLimit)
    : dpcp(dpcpRef),
      stab(),
      exactTimeLimit(exactTimeLimit),
      cxenv(),
      cxmodel(cxenv),
      y(cxenv, num_vertices(dpcp.get_graph())),
      w(cxenv, dpcp.get_nQ()),
      cxcons(cxenv),
      cplex(cxenv),
      cb(dpcpRef, stab, y, w),
      contextMask(0),
      mwis_env(NULL),
      mwis_pi(NULL),
      ecount(0),
      elist(NULL) {
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
  for (auto v : boost::make_iterator_range(vertices(dpcp.get_graph()))) {
    char name[100];
    snprintf(name, sizeof(name), "y_%ld_%ld", dpcp.get_P_part(v),
             dpcp.get_Q_part(v));
    y[dpcp.get_current_id(v)] = IloBoolVar(cxenv, name);
  }
  for (size_t qj = 0; qj < dpcp.get_nQ(); ++qj) {
    char name[100];
    snprintf(name, sizeof(name), "w_%ld", qj);
    w[qj] = IloBoolVar(cxenv, name);
  }

  // Objective
  IloExpr obj(cxenv);
  for (auto v : boost::make_iterator_range(vertices(dpcp.get_graph())))
    obj += y[dpcp.get_current_id(v)] * 1.0;
  for (size_t iB = 0; iB < dpcp.get_nQ(); ++iB) obj += w[iB] * (-1.0);
  cxobj = IloMaximize(cxenv, obj);
  cxmodel.add(cxobj);
  obj.end();

  // Constraints
  // (1) \sum_{j in Q: (i,j) in V} y_{i,j} <= 1, for all i in P
  for (size_t idP = 0; idP < dpcp.get_nP(); ++idP) {
    IloExpr restr(cxenv);
    for (Vertex v : dpcp.get_P()[idP]) restr += y[dpcp.get_current_id(v)];
    cxcons.add(restr <= 1);
    restr.end();
  }

  // (2) y_{i,j} + y_{i',j} <= w_j, for all ((i,j),(i',j)) in E such that i !=
  // i'
  for (auto e : boost::make_iterator_range(edges(dpcp.get_graph()))) {
    auto u = source(e, dpcp.get_graph());
    auto v = target(e, dpcp.get_graph());
    if (dpcp.get_P_part(u) == dpcp.get_P_part(v)) continue;
    if (dpcp.get_Q_part(u) != dpcp.get_Q_part(v)) continue;
    IloExpr restr(cxenv);
    restr += y[dpcp.get_current_id(u)] + y[dpcp.get_current_id(v)] -
             w[dpcp.get_Q_part(u)];
    cxcons.add(restr <= 0);
    restr.end();
  }

  // (3) y_{i,j} + \sum_{j' != j: (i',j') in N(i,j)} y_{i',j'} <= 1, for all
  // (i,j) in V, i != i'
  for (auto v : boost::make_iterator_range(vertices(dpcp.get_graph()))) {
    size_t pi = dpcp.get_P_part(v);
    size_t qj = dpcp.get_Q_part(v);
    for (size_t idP2 = 0; idP2 < dpcp.get_nP(); ++idP2) {
      if (idP2 == pi) continue;
      std::list<Vertex> neighbors;
      for (auto v2 : dpcp.get_P()[idP2])
        if (dpcp.get_Q_part(v2) != qj && edge(v, v2, dpcp.get_graph()).second)
          neighbors.push_back(v2);
      if (!neighbors.empty()) {
        IloExpr restr(cxenv);
        restr += y[dpcp.get_current_id(v)];
        for (auto v2 : neighbors) {
          restr += y[dpcp.get_current_id(v2)];
        }
        cxcons.add(restr <= 1);
        restr.end();
      }
    }
  }

  // (4) y_{i,j} <= w_j, for all (i,j) in V
  for (auto v : boost::make_iterator_range(vertices(dpcp.get_graph()))) {
    IloExpr restr(cxenv);
    restr += y[dpcp.get_current_id(v)] - w[dpcp.get_Q_part(v)];
    cxcons.add(restr <= 0);
    restr.end();
  }

  cxmodel.add(cxcons);

  // Re-export model
  cplex.extract(cxmodel);

  // Set parameters
  cplex.setDefaults();
  cplex.setParam(IloCplex::Param::Parallel, 1);  // Deterministic mode
  cplex.setParam(IloCplex::Param::Threads, 1);   // Single thread
  cplex.setOut(cxenv.getNullStream());
  cplex.setWarning(cxenv.getNullStream());

  // Now we get to setting up the generic callback.
  // We instantiate a ThresholdCallback and set the contextMask parameter.
  contextMask |= IloCplex::Callback::Context::Id::Candidate;

  // If contextMask is not zero we add the callback.
  if (contextMask != 0) cplex.use(&cb, contextMask);
}

void PricingEnv::mwis_init() {
  // Initalize stable environment
  COLORstable_initenv(&mwis_env, NULL, 0);

  // Intialize vectors of weights
  mwis_pi =
      (COLORNWT*)COLOR_SAFE_MALLOC(num_vertices(dpcp.get_graph()), COLORNWT);

  // Initialize edge array
  elist = (int*)malloc(sizeof(int) * 2 * num_edges(dpcp.get_graph()));
  for (auto e : boost::make_iterator_range(edges(dpcp.get_graph()))) {
    elist[2 * ecount] = dpcp.get_current_id(source(e, dpcp.get_graph()));
    elist[2 * ecount++ + 1] = dpcp.get_current_id(target(e, dpcp.get_graph()));
  }
}

/* Adaptation of  COLOR_double2COLORNWT from exactcolors */
int double2COLORNWT(COLORNWT nweights[], COLORNWT* scalef,
                    const std::vector<double>& dbl_weights) {
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
    assert(weight <= (double)COLORNWT_MAX);
    nweights[i] = (COLORNWT)weight;
  }
  return 0;
}

std::list<std::pair<StableEnv, PRICING_STATE>> PricingEnv::mwis_P_Q_solve(
    IloNumArray& dualsP, IloNumArray& dualsQ) {
  // Initialize local variables
  COLORset* newsets = NULL;
  int nnewsets = 0;
  COLORNWT mwis_pi_scalef = 1;
  MWISenv* mwis_env2 = NULL;
  COLORNWT* mwis_pi2 = NULL;
  COLORstable_initenv(&mwis_env2, NULL, 0);

  // Subgraph induced by vertices with positive weight

  // vertex mapping for graph to subgraph
  std::vector<int> vmap(num_vertices(dpcp.get_graph()), -1);
  // vertex mapping from subgraph to graph
  std::vector<int> invmap;
  std::vector<double> weights2;
  size_t i = 0;
  for (auto v : boost::make_iterator_range(vertices(dpcp.get_graph()))) {
    size_t pi = dpcp.get_P_part(v);
    size_t qj = dpcp.get_Q_part(v);
    double w = dualsP[pi] - dualsQ[qj];
    if (w > PRICING_EPSILON) {
      vmap[dpcp.get_current_id(v)] = i;
      invmap.push_back(dpcp.get_current_id(v));
      weights2.push_back(w);
      ++i;
    }
  }

  // Recompute edge list according to the previous indices
  // Ignore edges that have as endpoint a vertex with a negative weight
  int ecount2 = 0;
  int* elist2 = (int*)malloc(sizeof(int) * 2 * num_edges(dpcp.get_graph()));
  for (auto e : boost::make_iterator_range(edges(dpcp.get_graph()))) {
    auto v = source(e, dpcp.get_graph());
    auto u = target(e, dpcp.get_graph());
    if (vmap[dpcp.get_current_id(v)] == -1 ||
        vmap[dpcp.get_current_id(u)] == -1)
      continue;
    elist2[2 * ecount2] = vmap[dpcp.get_current_id(v)];
    elist2[2 * ecount2++ + 1] = vmap[dpcp.get_current_id(u)];
  }

  // Scale weights
  mwis_pi2 = (COLORNWT*)COLOR_SAFE_MALLOC(weights2.size(), COLORNWT);
  double2COLORNWT(mwis_pi2, &mwis_pi_scalef, weights2);

  if (ecount2 == 0) {
    // Edge-less graphs raise error in COLORstable_wrapper
    // So, they are manually solved
    nnewsets = 1;
    newsets = (COLORset*)malloc(sizeof(COLORset) * nnewsets);
    newsets[0].next = NULL;
    newsets[0].count = weights2.size();
    newsets[0].members = (int*)malloc(sizeof(int) * newsets[0].count);
    for (int i = 0; i < newsets[0].count; ++i) newsets[0].members[i] = i;
  } else
    // Solve the MWIS problem
    COLORstable_wrapper(&mwis_env2, &newsets, &nnewsets, weights2.size(),
                        ecount2, elist2, mwis_pi2, mwis_pi_scalef, 0, 0, 2);

  assert(nnewsets > 0);

  std::list<std::pair<StableEnv, PRICING_STATE>> ret;
  for (int set_i = 0; set_i < nnewsets; ++set_i) {
    // Recover stable set
    StableEnv st;
    for (int j = 0; j < newsets[set_i].count; ++j) {
      int v = invmap[newsets[set_i].members[j]];
      Vertex vv = vertex(v, dpcp.get_graph());
      st.stable.push_back(vv);
      size_t pi = dpcp.get_P_part(vv);
      size_t qj = dpcp.get_Q_part(vv);
      if (st.ps.insert(pi).second) st.cost += dualsP[pi];
      if (st.qs.insert(qj).second) st.cost -= dualsQ[qj];
    }

    // Maximalize stable set
    for (size_t qj : st.qs)
      for (Vertex v : dpcp.get_Q()[qj]) {
        // Check adjacencies
        bool ok = true;
        for (Vertex u : st.stable)
          if (u == v || edge(u, v, dpcp.get_graph()).second) ok = false;
        if (!ok) continue;
        size_t pi = dpcp.get_P_part(v);
        if (st.ps.insert(pi).second) {
          st.stable.push_back(v);
          st.cost += dualsP[pi];
        }
      }

    // Classify stable set
    PRICING_STATE state;
    if (st.cost > 1 + PRICING_EPSILON)
      state = PRICING_STABLE_FOUND;
    else
      state = PRICING_STABLE_NOT_FOUND;

    ret.push_back(std::make_pair(st, state));
  }

  // Free memory
  for (int i = 0; i < nnewsets; ++i) free(newsets[i].members);
  free(newsets);
  free(elist2);
  free(mwis_pi2);
  COLORstable_freeenv(&mwis_env2);

  return ret;
}

std::pair<StableEnv, PRICING_STATE> PricingEnv::mwis_P_solve(
    IloNumArray& dualsP, IloNumArray& dualsQ) {
  // Initialize local variables
  COLORset* newsets = NULL;
  int nnewsets = 0;
  COLORNWT mwis_pi_scalef = 1;

  // Reset stable
  stab.stable.clear();
  stab.ps.clear();
  stab.qs.clear();
  stab.cost = 0.0;

  // Compute vertex weight array
  std::vector<double> weights(num_vertices(dpcp.get_graph()));
  for (auto v : boost::make_iterator_range(vertices(dpcp.get_graph()))) {
    double w = dualsP[dpcp.get_P_part(v)];
    assert(w >= -PRICING_EPSILON);
    if (w < PRICING_EPSILON)
      weights[dpcp.get_current_id(v)] = 0.0;
    else
      weights[dpcp.get_current_id(v)] = w;
  }

  // Scale weights
  double2COLORNWT(mwis_pi, &mwis_pi_scalef, weights);
  mwis_pi_scalef = INT_MAX;  // Force optimality

  // Assert for non-negative costs
  for (size_t i = 0; i < weights.size(); ++i) {
    if (weights[i] < 0) {
      std::cerr << "Error: negative weight in mwis_P_solve: w[" << i
                << "] = " << weights[i] << std::endl;
      exit(1);
    }
  }

  // Solve the MWIS problem up to optimality
  COLORstable_wrapper(&mwis_env, &newsets, &nnewsets,
                      num_vertices(dpcp.get_graph()), ecount, elist, mwis_pi,
                      mwis_pi_scalef, 0, 0, 2);

  assert(nnewsets > 0);

  // Recover stable set
  // newsets[0] has an optimal stable set
  double w = 0.0;  // original weight
  for (int j = 0; j < newsets[0].count; ++j) {
    int v = newsets[0].members[j];
    w += weights[v];
    Vertex vv = vertex(v, dpcp.get_graph());
    stab.stable.push_back(vv);
    size_t pi = dpcp.get_P_part(vv);
    size_t qj = dpcp.get_Q_part(vv);
    if (stab.ps.insert(pi).second) stab.cost += dualsP[pi];
    if (stab.qs.insert(qj).second) stab.cost -= dualsQ[qj];
  }

  // Free memory
  for (int i = 0; i < nnewsets; ++i) free(newsets[i].members);
  free(newsets);

  // Classify stable set
  PRICING_STATE state;
  if (w <= 1 + PRICING_EPSILON)
    state = PRICING_STABLE_NOT_EXIST;
  else if (stab.cost <= 1 + PRICING_EPSILON)
    state = PRICING_STABLE_NOT_FOUND;
  else
    state = PRICING_STABLE_FOUND;

  return std::make_pair(stab, state);
}

std::pair<StableEnv, PRICING_STATE> PricingEnv::exact_solve(
    IloNumArray& dualsP, IloNumArray& dualsQ) {
  // Update objective coefficients
  IloNumArray y_coefs(cxenv, num_vertices(dpcp.get_graph()));
  for (auto v : boost::make_iterator_range(vertices(dpcp.get_graph())))
    y_coefs[dpcp.get_current_id(v)] = dualsP[dpcp.get_P_part(v)];
  IloNumArray w_coefs(cxenv, dpcp.get_nQ());
  for (size_t qj = 0; qj < dpcp.get_nQ(); ++qj) w_coefs[qj] = -dualsQ[qj];
  cxobj.setLinearCoefs(y, y_coefs);
  cxobj.setLinearCoefs(w, w_coefs);
  cxmodel.add(cxobj);

  // Reset stable
  stab.stable.clear();
  stab.ps.clear();
  stab.qs.clear();
  stab.cost = 0.0;

  // Solve
  cplex.extract(cxmodel);
  cplex.setParam(IloCplex::Param::TimeLimit, exactTimeLimit);
  cplex.solve();

  // Get final state
  PRICING_STATE state;
  switch (cplex.getCplexStatus()) {
    case IloCplex::CplexStatus::OptimalTol:
    case IloCplex::CplexStatus::Optimal: {
      std::cout << "Exact pricing: optimal solution found with objective value "
                << cplex.getObjValue() << std::endl;
      // The optimal value is <= THRESHOLD = 1.1
      // Recover optimal solution
      IloNumArray valY(cxenv, num_vertices(dpcp.get_graph()));
      cplex.getValues(y, valY);
      IloNumArray valW(cxenv, dpcp.get_nQ());
      cplex.getValues(w, valW);
      for (size_t qj = 0; qj < dpcp.get_nQ(); ++qj) {
        if (valW[qj] < 0.5) continue;
        // Postprocessing: force w_j <= \sum_{v in Q[j]} y_v
        size_t sum_y = 0;
        for (auto v : dpcp.get_Q()[qj])
          if (valY[dpcp.get_current_id(v)] > 0.5) sum_y += 1;
        if (sum_y == 0)
          valW[qj] = 0.0;
        else
          stab.qs.insert(qj);
      }
      for (auto v : boost::make_iterator_range(vertices(dpcp.get_graph())))
        if (valY[dpcp.get_current_id(v)] > 0.5) {
          stab.stable.push_back(v);
          stab.ps.insert(dpcp.get_P_part(v));
        }
      stab.cost = cplex.getBestObjValue();
      // Classify optimal solution
      if (stab.cost > 1 + PRICING_EPSILON)
        state = PRICING_STABLE_FOUND;
      else
        state = PRICING_STABLE_NOT_EXIST;
    } break;
    case IloCplex::CplexStatus::AbortUser:
      // Exit with abortion means that the current value is > THRESHOLD = 1.1
      state = PRICING_STABLE_FOUND;
      // Try to maximalize the stable set
      for (Vertex v : boost::make_iterator_range(vertices(dpcp.get_graph()))) {
        size_t pi = dpcp.get_P_part(v);
        size_t qj = dpcp.get_Q_part(v);
        if (dualsP[pi] < PRICING_EPSILON || stab.ps.contains(pi)) continue;
        if (!stab.qs.contains(qj)) continue;
        auto it = std::find_if(
            stab.stable.begin(), stab.stable.end(),
            [v, this](auto w) { return edge(v, w, dpcp.get_graph()).second; });
        if (it != stab.stable.end()) continue;
        stab.stable.push_back(v);
        stab.ps.insert(pi);
        stab.cost += dualsP[pi];
      }
      std::cout << "Exact pricing: threshold exceeded, with value " << stab.cost
                << std::endl;
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

std::pair<StableEnv, PRICING_STATE> PricingEnv::heur_solve(IloNumArray& dualsP,
                                                           IloNumArray& dualsQ,
                                                           double alpha) {
  // Reset stable
  stab.stable.clear();
  stab.ps.clear();
  stab.qs.clear();
  stab.cost = 0.0;

  // Candidates
  std::list<std::pair<double, Vertex>> candidates;
  double max_cost = std::numeric_limits<double>::lowest();
  double min_cost = std::numeric_limits<double>::max();
  for (auto v : boost::make_iterator_range(vertices(dpcp.get_graph()))) {
    size_t pi = dpcp.get_P_part(v);
    size_t qj = dpcp.get_Q_part(v);
    double cost_p = dualsP[pi];
    double cost_q = dualsQ[qj];
    double cost = cost_p - cost_q;
    candidates.push_back(std::make_pair(cost, v));
    if (cost > max_cost) max_cost = cost;
    if (cost < min_cost) min_cost = cost;
  }

  while (!candidates.empty()) {
    // Build RCL list
    std::list<std::pair<double, Vertex>> rcl;
    for (auto& [cost, v] : candidates)
      if (cost >= max_cost - alpha * (max_cost - min_cost) - PRICING_EPSILON)
        rcl.push_back(std::make_pair(cost, v));

    // Choose a random candidate from the RCL
    size_t random_index = rand_int(rng) % rcl.size();
    auto it_v = rcl.begin();
    std::advance(it_v, random_index);
    Vertex v = it_v->second;
    size_t pi = dpcp.get_P_part(v);
    size_t qj = dpcp.get_Q_part(v);

    // Add the best candidate to the stable set
    stab.stable.push_back(v);
    stab.ps.insert(pi);
    stab.qs.insert(qj);
    stab.cost += it_v->first;

    // Remove and update weights in the candidate list
    max_cost = std::numeric_limits<double>::lowest();
    min_cost = std::numeric_limits<double>::max();
    for (auto it = candidates.begin(); it != candidates.end();) {
      Vertex u = it->second;
      size_t pu = dpcp.get_P_part(u);
      size_t qu = dpcp.get_Q_part(u);
      if (pu == pi || edge(u, v, dpcp.get_graph()).second) {
        it = candidates.erase(it);
        continue;
      }
      if (qu == qj) it->first = dualsP[pu];
      if (it->first > max_cost) max_cost = it->first;
      if (it->first < min_cost) min_cost = it->first;
      ++it;
    }
  }

  PRICING_STATE state;
  if (stab.cost > 1 + PRICING_EPSILON)
    state = PRICING_STABLE_FOUND;
  else
    state = PRICING_STABLE_NOT_FOUND;

  return std::make_pair(stab, state);
}
