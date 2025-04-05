#include "pricing.hpp"
#include <limits>

std::uniform_int_distribution<rng_type::result_type>
    udist(0, std::numeric_limits<rng_type::result_type>::max());
std::uniform_real_distribution<double> cdist(0.0, 1.0);

void ThresholdCallback::check_thresolhd(
    const IloCplex::Callback::Context &context) {
  if (context.isCandidatePoint())
    if (context.getCandidateObjective() > THRESHOLD) {
      // Save stable set
      IloNumArray valY(context.getEnv(), num_vertices(in.graph));
      context.getCandidatePoint(y, valY);
      IloNumArray valWa(context.getEnv(), in.nA);
      context.getCandidatePoint(wa, valWa);
      IloNumArray valWb(context.getEnv(), in.nB);
      context.getCandidatePoint(wb, valWb);
      for (auto v : boost::make_iterator_range(vertices(in.graph)))
        if (valY[v] > 0.5)
          stab.stable.push_back(v);
      for (size_t iA = 0; iA < in.nA; ++iA)
        if (valWa[iA] > 0.5)
          stab.as.push_back(in.idA2TyA[iA]);
      for (size_t iB = 0; iB < in.nB; ++iB)
        if (valWb[iB] > 0.5)
          stab.bs.push_back(in.idB2TyB[iB]);
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
      wa(cxenv, in.nA), wb(cxenv, in.nB), cxcons(cxenv), cplex(cxenv),
      cb(in, stab, y, wa, wb), contextMask(0), rng() {
  exact_init();
  // seed rng:
  rng_type::result_type const seedval = 0;
  rng.seed(seedval);
}

PricingEnv::~PricingEnv() {
  // End CPLEX variables
  cplex.end();
  cxcons.end();
  y.end();
  wa.end();
  wb.end();
  cxobj.end();
  cxmodel.end();
  cxenv.end();
}

void PricingEnv::exact_init() {
  // Variables
  for (auto v : boost::make_iterator_range(vertices(in.graph))) {
    char name[100];
    snprintf(name, sizeof(name), "y_%ld_%ld", in.graph[v].first,
             in.graph[v].second);
    y[v] = IloBoolVar(cxenv, name);
  }
  for (size_t iA = 0; iA < in.nA; ++iA) {
    char name[100];
    snprintf(name, sizeof(name), "wa_%ld", iA);
    wa[iA] = IloBoolVar(cxenv, name);
  }
  for (size_t iB = 0; iB < in.nB; ++iB) {
    char name[100];
    snprintf(name, sizeof(name), "wb_%ld", iB);
    wb[iB] = IloBoolVar(cxenv, name);
  }

  // Objective
  IloExpr obj(cxenv);
  for (size_t iA = 0; iA < in.nA; ++iA)
    obj += wa[iA] * 1.0;
  for (size_t iB = 0; iB < in.nB; ++iB)
    obj += wb[iB] * (-1.0);
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

  // (3a) y_a_b <= wa_a, for all (a,b) \in V
  for (auto v : boost::make_iterator_range(vertices(in.graph))) {
    IloExpr restr(cxenv);
    restr += y[v] - wa[in.tyA2idA[in.graph[v].first]];
    cxcons.add(restr <= 0);
    restr.end();
  }
  for (size_t iA = 0; iA < in.nA; ++iA) {
    IloExpr restr(cxenv);
    restr += wa[iA];
    for (auto v : in.snd[iA])
      restr -= y[v];
    cxcons.add(restr <= 0);
    restr.end();
  }

  // (3b) y_a_b <= wb_b, for all (a,b) \in V
  for (auto v : boost::make_iterator_range(vertices(in.graph))) {
    IloExpr restr(cxenv);
    restr += y[v] - wb[in.tyB2idB[in.graph[v].second]];
    cxcons.add(restr <= 0);
    restr.end();
  }
  for (size_t iB = 0; iB < in.nB; ++iB) {
    IloExpr restr(cxenv);
    restr += wb[iB];
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
  cplex.setParam(IloCplex::Param::TimeLimit, TIMELIMIT);
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

std::pair<StableEnv, PRICING_STATE>
PricingEnv::exact_solve(IloNumArray &dualsA, IloNumArray &dualsB) {
  // Update objective coefficients
  cxobj.setLinearCoefs(wa, dualsA);
  cxobj.setLinearCoefs(wb, dualsB);
  cxmodel.add(cxobj);

  // Reset stable
  stab.stable.clear();
  stab.as.clear();
  stab.bs.clear();

  // Solve
  cplex.extract(cxmodel);
  cplex.solve();

  // Get final state
  PRICING_STATE state;
  switch (cplex.getCplexStatus()) {
  case IloCplex::CplexStatus::Optimal:
    state = PRICING_OPTIMAL;
    break;
  case IloCplex::CplexStatus::AbortUser:
    state = PRICING_FEASIBLE;
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

  // Recover solution
  if (state == PRICING_OPTIMAL) {
    IloNumArray valY(cxenv, num_vertices(in.graph));
    cplex.getValues(y, valY);
    IloNumArray valWa(cxenv, in.nA);
    cplex.getValues(wa, valWa);
    IloNumArray valWb(cxenv, in.nB);
    cplex.getValues(wb, valWb);
    for (auto v : boost::make_iterator_range(vertices(in.graph)))
      if (valY[v] > 0.5)
        stab.stable.push_back(v);
    for (size_t iA = 0; iA < in.nA; ++iA)
      if (valWa[iA] > 0.5)
        stab.as.push_back(in.idA2TyA[iA]);
    for (size_t iB = 0; iB < in.nB; ++iB)
      if (valWb[iB] > 0.5)
        stab.bs.push_back(in.idB2TyB[iB]);
    stab.cost = cplex.getBestObjValue();
  }

  return std::make_pair(stab, state);
}

void PricingEnv::heur_init(IloNumArray &dualsA, IloNumArray &dualsB) {

  // Create a list with (delta_a - mu_b, v, a, best_b) for all a in A
  heurCandidates.clear();
  for (size_t iA = 0; iA < in.nA; ++iA) {
    // Find best b
    assert(!in.snd[iA].empty());
    size_t best_v = 0;
    TypeB best_b = 0;
    double best_weight = -std::numeric_limits<double>::max();
    for (auto v : in.snd[iA]) {
      size_t b = in.graph[v].second;
      double weight = dualsA[iA] + dualsB[in.tyB2idB[b]];
      if (weight > best_weight) {
        best_weight = weight;
        best_b = b;
        best_v = v;
      }
    }
    heurCandidates.push_back(
        std::make_tuple(best_weight, best_v, in.idA2TyA[iA], best_b));
  }
  assert(!heurCandidates.empty());

  // Sort candidates
  heurCandidates.sort(std::greater<>());
}

void PricingEnv::try_stable_add(double weight, size_t v, TypeA a, TypeB b,
                                IloNumArray &dualsA, std::vector<bool> &used,
                                std::set<TypeB> &bs) {
  if (used[in.tyA2idA[a]])
    return;

  if (std::find_if(stab.stable.begin(), stab.stable.end(), [v, this](auto w) {
        return edge(v, w, in.graph).second;
      }) != stab.stable.end())
    return;

  // Add v = (a,b) to the stable set
  stab.stable.push_back(v);
  stab.as.push_back(a);
  stab.cost += weight;
  bs.insert(b);
  used[in.tyA2idA[a]] = true;

  // Use b to represent other vertices
  std::shuffle(in.fst[in.tyB2idB[b]].begin(), in.fst[in.tyB2idB[b]].end(), rng);
  for (size_t u : in.fst[in.tyB2idB[b]]) {

    if (cdist(rng) < 0.5)
      continue;

    if (edge(v, u, in.graph).second)
      continue;

    TypeA au = in.graph[u].first;
    if (used[in.tyA2idA[au]])
      continue;

    if (std::find_if(stab.stable.begin(), stab.stable.end(), [u, this](auto w) {
          return edge(u, w, in.graph).second;
        }) != stab.stable.end())
      continue;

    // Add u = (au,b) to the stable set
    stab.stable.push_back(u);
    stab.as.push_back(au);
    stab.cost += dualsA[in.tyA2idA[au]];
    used[in.tyA2idA[au]] = true;
  }
}

std::pair<StableEnv, PRICING_STATE>
PricingEnv::heur_solve(IloNumArray &dualsA, IloNumArray &dualsB, TypeA first) {

  // Reset stable
  stab.stable.clear();
  stab.as.clear();
  stab.bs.clear();
  stab.cost = 0.0;

  std::vector<bool> used(in.nA, false);
  std::set<TypeB> bs;

  // Add the first vertex to start the sable set
  auto it = std::find_if(heurCandidates.begin(), heurCandidates.end(),
                         [first](auto p) { return std::get<2>(p) == first; });
  assert(it != heurCandidates.end());
  try_stable_add(std::get<0>(*it), std::get<1>(*it), std::get<2>(*it),
                 std::get<3>(*it), dualsA, used, bs);

  for (auto [weight, v, a, b] : heurCandidates) {
    // Try to add candidate to the stable set
    try_stable_add(weight, v, a, b, dualsA, used, bs);
  }

  // Update bs
  for (size_t b : bs)
    stab.bs.push_back(b);

  return std::make_pair(stab, PRICING_FEASIBLE);
}