#include "pricing.hpp"

#define TIMELIMIT 300.0 // = 5 minutes

std::tuple<StableEnv, PRICING_STATE> exact_solve(GraphEnv &in,
                                                 IloNumArray &duals) {

  // CPLEX enviroment variables
  IloEnv cxenv;
  IloModel cxmodel(cxenv);

  // Variables
  IloNumVarArray y(cxenv, num_vertices(in.graph));
  IloNumVarArray w(cxenv, in.nB);
  for (auto v : boost::make_iterator_range(vertices(in.graph))) {
    char name[100];
    snprintf(name, sizeof(name), "y_%ld_%ld", in.graph[v].first,
             in.graph[v].second);
    y[v] = IloBoolVar(cxenv, name);
  }
  for (size_t ib = 0; ib < in.nB; ++ib) {
    char name[100];
    snprintf(name, sizeof(name), "w_%ld", ib);
    w[ib] = IloBoolVar(cxenv, name);
  }

  // Objective
  IloExpr fobj(cxenv, 0);
  for (auto v : boost::make_iterator_range(vertices(in.graph)))
    fobj += y[v] * duals[in.tyA2idA[in.graph[v].first]];
  for (size_t iB = 0; iB < in.nB; ++iB)
    fobj -= w[iB] * duals[in.nA + iB];
  cxmodel.add(IloMaximize(cxenv, fobj));

  // Constraints
  IloConstraintArray cxcons(cxenv);

  // (1) \sum_{b \in B: (a,b) \in V} y_a_b <= 1, for all a \in A
  for (size_t ia = 0; ia < in.nA; ++ia) {
    IloExpr restr(cxenv);
    for (size_t v : in.snd[ia])
      restr += y[v];
    cxcons.add(restr <= 1);
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
  }

  // (3) y_a_b <= w_b, for all (a,b) \in V
  for (auto v : boost::make_iterator_range(vertices(in.graph))) {
    IloExpr restr(cxenv);
    restr += y[v] - w[in.tyB2idB[in.graph[v].second]];
    cxcons.add(restr <= 0);
  }

  cxmodel.add(cxcons);

  // CPLEX
  IloCplex cplex(cxmodel);

  // Set parameters
  cplex.setDefaults();
  cplex.setParam(IloCplex::Param::TimeLimit, TIMELIMIT);
  cplex.setParam(IloCplex::Param::Parallel, 1); // Deterministic mode
  cplex.setParam(IloCplex::Param::Threads, 1);  // Single thread
  cplex.setOut(cxenv.getNullStream());
  cplex.setWarning(cxenv.getNullStream());

  // Solve
  cplex.solve();

  // Get final state
  PRICING_STATE state;
  switch (cplex.getCplexStatus()) {
  case IloCplex::CplexStatus::Optimal:
    state = PRICING_OPTIMAL;
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
  StableEnv stab;
  if (state == PRICING_OPTIMAL) {
    for (auto v : boost::make_iterator_range(vertices(in.graph)))
      if (cplex.getValue(y[v]) > 0.5) {
        stab.stable.push_back(v);
        TypeA a = in.graph[v].first;
        stab.cost += duals[in.tyA2idA[a]];
      }
    for (size_t iB = 0; iB < in.nB; ++iB)
      if (cplex.getValue(w[iB]) > 0.5) {
        stab.bs.push_back(in.idB2TyB[iB]);
        stab.cost -= duals[in.nA + iB];
      }
  }

  // End CPLEX variables
  cplex.end();
  cxcons.end();
  y.end();
  w.end();
  fobj.end();
  cxmodel.end();
  cxenv.end();

  return std::make_tuple(stab, state);
}