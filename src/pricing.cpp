#include "pricing.hpp"
#include "random.hpp"
#include <cfloat>
#include <limits>

// Generic callback: abort when the threshold is exceeded.
// Implementation of early stop.
void ThresholdCallback::check_threshold(
    const IloCplex::Callback::Context &context) {
  if (context.isCandidatePoint())
    if (context.getCandidateObjective() > THRESHOLD) {
      // Save stable set
      IloNumArray valY(context.getEnv(), num_vertices(in.graph));
      context.getCandidatePoint(y, valY);
      IloNumArray valWb(context.getEnv(), in.nB);
      context.getCandidatePoint(w, valWb);
      for (size_t iB = 0; iB < in.nB; ++iB) {
        if (valWb[iB] < 0.5)
          continue;
        TypeB b = in.idB2TyB[iB];
        // Postprocessing: force w_b \leq \sum_{v \in B^b} y_v
        size_t sum_y = 0;
        for (auto v : in.Vb[b])
          if (valY[in.getId[v]] > 0.5)
            sum_y += 1;
        if (sum_y == 0)
          valWb[iB] = 0.0;
        else
          stab.bs.insert(in.idB2TyB[iB]);
      }
      for (auto v : boost::make_iterator_range(vertices(in.graph)))
        if (valY[in.getId[v]] > 0.5) {
          stab.stable.push_back(v);
          stab.as.insert(in.graph[v].first);
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
void ThresholdCallback::invoke(const IloCplex::Callback::Context &context) {
  if (context.inCandidate())
    check_threshold(context);
}

PricingEnv::PricingEnv(GraphEnv &in, double exactTimeLimit)
    : in(in), stab(), exactTimeLimit(exactTimeLimit), cxenv(), cxmodel(cxenv),
      y(cxenv, num_vertices(in.graph)), w(cxenv, in.nB), cxcons(cxenv),
      cplex(cxenv), cb(in, stab, y, w), contextMask(0), mwis_env(NULL),
      mwis_pi(NULL), ecount(0), elist(NULL) {
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
    y[in.getId[v]] = IloBoolVar(cxenv, name);
  }
  for (size_t iB = 0; iB < in.nB; ++iB) {
    char name[100];
    snprintf(name, sizeof(name), "w_%ld", iB);
    w[iB] = IloBoolVar(cxenv, name);
  }

  // Objective
  IloExpr obj(cxenv);
  for (auto v : boost::make_iterator_range(vertices(in.graph)))
    obj += y[in.getId[v]] * 1.0;
  for (size_t iB = 0; iB < in.nB; ++iB)
    obj += w[iB] * (-1.0);
  cxobj = IloMaximize(cxenv, obj);
  cxmodel.add(cxobj);
  obj.end();

  // Constraints
  // (1) \sum_{b \in B: (a,b) \in V} y_a_b <= 1, for all a \in A
  for (size_t ia = 0; ia < in.nA; ++ia) {
    IloExpr restr(cxenv);
    for (Vertex v : in.snd[ia])
      restr += y[in.getId[v]];
    cxcons.add(restr <= 1);
    restr.end();
  }

  // (2) y_a_b + y_a'_b <= w_b, for all ((a,b),(a',b)) \in E such that a != a'
  for (auto e : boost::make_iterator_range(edges(in.graph))) {
    auto u = source(e, in.graph);
    auto v = target(e, in.graph);
    if (in.graph[u].first == in.graph[v].first)
      continue;
    if (in.graph[u].second != in.graph[v].second)
      continue;
    IloExpr restr(cxenv);
    restr +=
        y[in.getId[u]] + y[in.getId[v]] - w[in.tyB2idB[in.graph[u].second]];
    cxcons.add(restr <= 0);
    restr.end();
  }

  // (3) y_a_b + \sum_{b' != b: (a',b') \in N(a,b)} y_a'_b' <= 1, for all
  // (a,b) \in V, a != a'
  for (auto v : boost::make_iterator_range(vertices(in.graph))) {
    TypeA a = in.graph[v].first;
    TypeB b = in.graph[v].second;
    for (size_t a2 = 0; a2 < in.nA; ++a2) {
      if (a2 == a)
        continue;
      std::list<Vertex> neighbors;
      for (auto v2 : in.snd[a2])
        if (in.graph[v2].second != b && edge(v, v2, in.graph).second)
          neighbors.push_back(v2);
      if (!neighbors.empty()) {
        IloExpr restr(cxenv);
        restr += y[in.getId[v]];
        for (auto v2 : neighbors) {
          restr += y[in.getId[v2]];
        }
        cxcons.add(restr <= 1);
        restr.end();
      }
    }
  }

  // (4) y_a_b <= w_b, for all (a,b) \in V
  for (auto v : boost::make_iterator_range(vertices(in.graph))) {
    IloExpr restr(cxenv);
    restr += y[in.getId[v]] - w[in.tyB2idB[in.graph[v].second]];
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

  // Initialize edge array
  elist = (int *)malloc(sizeof(int) * 2 * num_edges(in.graph));
  for (auto e : boost::make_iterator_range(edges(in.graph))) {
    elist[2 * ecount] = in.getId[source(e, in.graph)];
    elist[2 * ecount++ + 1] = in.getId[target(e, in.graph)];
  }
}

/* Adaptation of  COLOR_double2COLORNWT from exactcolors */
int double2COLORNWT(COLORNWT nweights[], COLORNWT *scalef,
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
    assert(weight <= (double)COLORNWT_MAX);
    nweights[i] = (COLORNWT)weight;
  }
  return 0;
}

std::list<std::pair<StableEnv, PRICING_STATE>>
PricingEnv::mwis_P_Q_solve(IloNumArray &dualsA, IloNumArray &dualsB) {

  // Initialize local variables
  COLORset *newsets = NULL;
  int nnewsets = 0;
  COLORNWT mwis_pi_scalef = 1;
  MWISenv *mwis_env2 = NULL;
  COLORNWT *mwis_pi2 = NULL;
  COLORstable_initenv(&mwis_env2, NULL, 0);

  // Subgraph induced by vertices with positive weight

  // vertex mapping for graph to subgraph
  std::vector<int> vmap(num_vertices(in.graph), -1);
  // vertex mapping from subgraph to graph
  std::vector<int> invmap;
  std::vector<double> weights2;
  size_t i = 0;
  for (auto v : boost::make_iterator_range(vertices(in.graph))) {
    auto [a, b, id] = in.graph[v];
    double w = dualsA[in.tyA2idA[a]] - dualsB[in.tyB2idB[b]];
    if (w > PRICING_EPSILON) {
      vmap[in.getId[v]] = i;
      invmap.push_back(in.getId[v]);
      weights2.push_back(w);
      ++i;
    }
  }

  // Recompute edge list according to the previous indices
  // Ignore edges that have as endpoint a vertex with a negative weight
  int ecount2 = 0;
  int *elist2 = (int *)malloc(sizeof(int) * 2 * num_edges(in.graph));
  for (auto e : boost::make_iterator_range(edges(in.graph))) {
    auto v = source(e, in.graph);
    auto u = target(e, in.graph);
    if (vmap[in.getId[v]] == -1 || vmap[in.getId[u]] == -1)
      continue;
    elist2[2 * ecount2] = vmap[in.getId[v]];
    elist2[2 * ecount2++ + 1] = vmap[in.getId[u]];
  }

  // Scale weights
  mwis_pi2 = (COLORNWT *)COLOR_SAFE_MALLOC(weights2.size(), COLORNWT);
  double2COLORNWT(mwis_pi2, &mwis_pi_scalef, weights2);

  if (ecount2 == 0) {
    // Edge-less graphs raise error in COLORstable_wrapper
    // So, they are manually solved
    nnewsets = 1;
    newsets = (COLORset *)malloc(sizeof(COLORset) * nnewsets);
    newsets[0].next = NULL;
    newsets[0].count = weights2.size();
    newsets[0].members = (int *)malloc(sizeof(int) * newsets[0].count);
    for (int i = 0; i < newsets[0].count; ++i)
      newsets[0].members[i] = i;
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
      Vertex vv = vertex(v, in.graph);
      st.stable.push_back(vv);
      auto [a, b, id] = in.graph[vv];
      if (st.as.insert(a).second)
        st.cost += dualsA[in.tyA2idA[a]];
      if (st.bs.insert(b).second)
        st.cost -= dualsB[in.tyB2idB[b]];
    }

    // Maximalize stable set
    for (TypeB b : st.bs)
      for (Vertex v : in.Vb[b]) {
        // Check adjacencies
        bool ok = true;
        for (Vertex u : st.stable)
          if (u == v || edge(u, v, in.graph).second)
            ok = false;
        if (!ok)
          continue;
        TypeA a = in.graph[v].first;
        if (st.as.insert(a).second) {
          st.stable.push_back(v);
          st.cost += dualsA[in.tyA2idA[a]];
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
  for (int i = 0; i < nnewsets; ++i)
    free(newsets[i].members);
  free(newsets);
  free(elist2);
  free(mwis_pi2);
  COLORstable_freeenv(&mwis_env2);

  return ret;
}

std::pair<StableEnv, PRICING_STATE>
PricingEnv::mwis_P_solve(IloNumArray &dualsA, IloNumArray &dualsB) {

  // Initialize local variables
  COLORset *newsets = NULL;
  int nnewsets = 0;
  COLORNWT mwis_pi_scalef = 1;

  // Reset stable
  stab.stable.clear();
  stab.as.clear();
  stab.bs.clear();
  stab.cost = 0.0;

  // Compute vertex weight array
  std::vector<double> weights(num_vertices(in.graph));
  for (auto v : boost::make_iterator_range(vertices(in.graph))) {
    double w = dualsA[in.tyA2idA[in.graph[v].first]];
    assert(w >= -PRICING_EPSILON);
    if (w < PRICING_EPSILON)
      weights[in.getId[v]] = 0.0;
    else
      weights[in.getId[v]] = w;
  }

  // Scale weights
  double2COLORNWT(mwis_pi, &mwis_pi_scalef, weights);
  mwis_pi_scalef = INT_MAX; // Force optimality

  // Assert for non-negative costs
  for (size_t i = 0; i < weights.size(); ++i) {
    if (weights[i] < 0) {
      std::cerr << "Error: negative weight in mwis_P_solve: w[" << i
                << "] = " << weights[i] << std::endl;
      exit(1);
    }
  }

  // Solve the MWIS problem up to optimality
  COLORstable_wrapper(&mwis_env, &newsets, &nnewsets, num_vertices(in.graph),
                      ecount, elist, mwis_pi, mwis_pi_scalef, 0, 0, 2);

  assert(nnewsets > 0);

  // Recover stable set
  // newsets[0] has an optimal stable set
  double w = 0.0; // original weight
  for (int j = 0; j < newsets[0].count; ++j) {
    int v = newsets[0].members[j];
    w += weights[v];
    Vertex vv = vertex(v, in.graph);
    stab.stable.push_back(vv);
    auto [a, b, id] = in.graph[vv];
    if (stab.as.insert(a).second)
      stab.cost += dualsA[in.tyA2idA[a]];
    if (stab.bs.insert(b).second)
      stab.cost -= dualsB[in.tyB2idB[b]];
  }

  // Free memory
  for (int i = 0; i < nnewsets; ++i)
    free(newsets[i].members);
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

std::pair<StableEnv, PRICING_STATE>
PricingEnv::exact_solve(IloNumArray &dualsA, IloNumArray &dualsB) {

  // Update objective coefficients
  IloNumArray y_coefs(cxenv, num_vertices(in.graph));
  for (auto v : boost::make_iterator_range(vertices(in.graph)))
    y_coefs[in.getId[v]] = dualsA[in.tyA2idA[in.graph[v].first]];
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
  cplex.setParam(IloCplex::Param::TimeLimit, exactTimeLimit);
  cplex.solve();

  // Get final state
  PRICING_STATE state;
  switch (cplex.getCplexStatus()) {
  case IloCplex::CplexStatus::OptimalTol:
  case IloCplex::CplexStatus::Optimal: {
    // The optimal value is <= THRESHOLD = 1.1
    // Recover optimal solution
    IloNumArray valY(cxenv, num_vertices(in.graph));
    cplex.getValues(y, valY);
    IloNumArray valW(cxenv, in.nB);
    cplex.getValues(w, valW);
    for (size_t iB = 0; iB < in.nB; ++iB) {
      if (valW[iB] < 0.5)
        continue;
      TypeB b = in.idB2TyB[iB];
      // Postprocessing: force w_b \leq \sum_{v \in B^b} y_v
      size_t sum_y = 0;
      for (auto v : in.Vb[b])
        if (valY[in.getId[v]] > 0.5)
          sum_y += 1;
      if (sum_y == 0)
        valW[iB] = 0.0;
      else
        stab.bs.insert(in.idB2TyB[iB]);
    }
    for (auto v : boost::make_iterator_range(vertices(in.graph)))
      if (valY[in.getId[v]] > 0.5) {
        stab.stable.push_back(v);
        stab.as.insert(in.graph[v].first);
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
    for (Vertex v : boost::make_iterator_range(vertices(in.graph))) {
      auto [a, b, id] = in.graph[v];
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

std::pair<StableEnv, PRICING_STATE>
PricingEnv::heur_solve(IloNumArray &dualsA, IloNumArray &dualsB, double alpha) {

  // Reset stable
  stab.stable.clear();
  stab.as.clear();
  stab.bs.clear();
  stab.cost = 0.0;

  // Candidates
  std::list<std::pair<double, Vertex>> candidates;
  double max_cost = std::numeric_limits<double>::lowest();
  double min_cost = std::numeric_limits<double>::max();
  for (auto v : boost::make_iterator_range(vertices(in.graph))) {
    auto [a, b, id] = in.graph[v];
    double cost_a = dualsA[in.tyA2idA[a]];
    double cost_b = dualsB[in.tyB2idB[b]];
    double cost = cost_a - cost_b;
    candidates.push_back(std::make_pair(cost, v));
    if (cost > max_cost)
      max_cost = cost;
    if (cost < min_cost)
      min_cost = cost;
  }

  while (!candidates.empty()) {

    // Build RCL list
    std::list<std::pair<double, Vertex>> rcl;
    for (auto &[cost, v] : candidates)
      if (cost >= max_cost - alpha * (max_cost - min_cost) - PRICING_EPSILON)
        rcl.push_back(std::make_pair(cost, v));

    // Choose a random candidate from the RCL
    size_t random_index = rand_int(rng) % rcl.size();
    auto it_v = rcl.begin();
    std::advance(it_v, random_index);
    Vertex v = it_v->second;
    auto [a, b, id] = in.graph[v];

    // Add the best candidate to the stable set
    stab.stable.push_back(v);
    stab.as.insert(a);
    stab.bs.insert(b);
    stab.cost += it_v->first;

    // Remove and update weights in the candidate list
    max_cost = std::numeric_limits<double>::lowest();
    min_cost = std::numeric_limits<double>::max();
    for (auto it = candidates.begin(); it != candidates.end();) {
      Vertex u = it->second;
      auto [au, bu, idu] = in.graph[u];
      if (au == a || edge(u, v, in.graph).second) {
        it = candidates.erase(it);
        continue;
      }
      if (bu == b)
        it->first = dualsA[in.tyA2idA[au]];
      if (it->first > max_cost)
        max_cost = it->first;
      if (it->first < min_cost)
        min_cost = it->first;
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