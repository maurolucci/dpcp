#ifndef __STATS_HPP__
#define __STATS_HPP__

#include <iostream>
#include <limits>
#include <string>

enum PRICING_STATE {
  PRICING_STABLE_FOUND,
  PRICING_STABLE_NOT_FOUND,
  PRICING_STABLE_NOT_EXIST,
  PRICING_READY,
  PRICING_TIME_EXCEEDED,
  PRICING_MEM_EXCEEDED,
  PRICING_OTHER,
};

enum LP_STATE {
  LP_UNSOLVED,
  LP_INFEASIBLE,
  LP_INTEGER,
  LP_FRACTIONAL,
  LP_TIME_EXCEEDED,
  LP_TIME_EXCEEDED_PR,
  LP_MEM_EXCEEDED,
  LP_MEM_EXCEEDED_PR,
  LP_INIT_FAIL,
};

enum STATE {
  OPTIMAL,
  FEASIBLE,
  INFEASIBLE,
  TIME_EXCEEDED,
  TIME_EXCEEDED_LP,
  TIME_EXCEEDED_PR,
  MEM_EXCEEDED,
  MEM_EXCEEDED_LP,
  MEM_EXCEEDED_PR,
  INIT_FAIL,
  UNKNOWN,
};

class Stats {

public:
  // Instance name
  std::string instance;
  // Solver name
  std::string solver;
  // Run number
  int run;
  // Number of vertices
  int nvertices;
  // Number of edges
  int nedges;
  // Cardinality of A
  int nA;
  // Cardinality of B
  int nB;
  // Number of variables
  int nvars;
  // Number of constraints
  int ncons;
  // Final state
  STATE state;
  // Total time
  double time;
  // Number of processed nodes
  int nodes;
  // Number of nodes left in the queue
  int nodesLeft;
  // Final lower bound
  double lb;
  // Final upper bound
  int ub;
  // Final optimality gap (in percentage)
  double gap;
  // Number of infeasible instances detected in the byp tree
  int ninfeasPrepro, ninfeasCheck, ninfeasAux;
  // Number of integer nodes detected in the byp tree
  int nint;
  // Number of GCP instances detected in the byp tree
  int ngcp;
  // Time spent on solving GCP instances
  double gcpTime;
  // Number of solutions found by heuristic and linear relaxation
  int nsolHeur, nsolLR;
  // Root lower bound
  double rootlb;
  // Root upper bound
  int rootub;
  // Time required for finding the initial solution at the root node
  double rootHeurTime;
  // Time required for feasibility check at the root node
  double rootFeasTime;
  // For each pricing method, number of columns added, number of calls, and
  // total time required at the root node
  int rootNColsPool, rootNColsHeur, rootNColsMwis1, rootNColsMwis2,
      rootNColsExact;
  int rootNCallsPool, rootNCallsHeur, rootNCallsMwis1, rootNCallsMwis2,
      rootNCallsExact;
  double rootTimePool, rootTimeHeur, rootTimeMwis1, rootTimeMwis2,
      rootTimeExact;
  // Time required for finding heuristic solutions in at other nodes
  double otherNodesHeurTime;
  // Number of calls to the feasibility check at other nodes
  int otherNodesFeasNCalls;
  // Time required for feasibility checks at other nodes
  double otherNodesFeasTime;
  // For each pricing method, number of columns added, number of calls, and
  // total time required at other nodes
  int otherNodesNColsPool, otherNodesNColsHeur, otherNodesNColsMwis1,
      otherNodesNColsMwis2, otherNodesNColsExact;
  int otherNodesNCallsPool, otherNodesNCallsHeur, otherNodesNCallsMwis1,
      otherNodesNCallsMwis2, otherNodesNCallsExact;
  double otherNodesTimePool, otherNodesTimeHeur, otherNodesTimeMwis1,
      otherNodesTimeMwis2, otherNodesTimeExact;

  Stats()
      : instance(""), solver(""), run(-1), nvertices(-1), nedges(-1), nA(-1),
        nB(-1), nvars(-1), ncons(-1), state(UNKNOWN), time(0.0), nodes(0),
        nodesLeft(0), lb(-1.0), ub(-1), gap(-1.0), ninfeasPrepro(0),
        ninfeasCheck(0), ninfeasAux(0), nint(0), ngcp(0), gcpTime(0.0),
        nsolHeur(0), nsolLR(0), rootlb(-1.0), rootub(-1), rootHeurTime(0.0),
        rootFeasTime(0.0), rootNColsPool(0), rootNColsHeur(0),
        rootNColsMwis1(0), rootNColsMwis2(0), rootNColsExact(0),
        rootNCallsPool(0), rootNCallsHeur(0), rootNCallsMwis1(0),
        rootNCallsMwis2(0), rootNCallsExact(0), rootTimePool(0.0),
        rootTimeHeur(0.0), rootTimeMwis1(0.0), rootTimeMwis2(0.0),
        rootTimeExact(0.0), otherNodesHeurTime(0.0), otherNodesFeasNCalls(0),
        otherNodesFeasTime(0.0), otherNodesNColsPool(0), otherNodesNColsHeur(0),
        otherNodesNColsMwis1(0), otherNodesNColsMwis2(0),
        otherNodesNColsExact(0), otherNodesNCallsPool(0),
        otherNodesNCallsHeur(0), otherNodesNCallsMwis1(0),
        otherNodesNCallsMwis2(0), otherNodesNCallsExact(0),
        otherNodesTimePool(0.0), otherNodesTimeHeur(0.0),
        otherNodesTimeMwis1(0.0), otherNodesTimeMwis2(0.0),
        otherNodesTimeExact(0.0) {}

  std::string get_state_as_str();
  void write_stats(std::ostream &file);
  void print_stats(std::ostream &file);
};

#endif // __STATS_HPP__