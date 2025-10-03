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
  std::string instance;
  std::string solver;
  int run;
  int nvertices;
  int nedges;
  int nA;
  int nB;
  int nvars;
  int ncons;
  STATE state;
  double time;
  int nodes;
  int nodesLeft;
  int initSol;
  double initSolTime;
  double rootval;
  double lb;
  int ub;
  double gap;
  int ninfeas;
  int ngcp;
  int poolSize;
  int nColsPool;
  int nColsHeur;
  int nColsMwis1;
  int nColsMwis2;
  int nColsExact;
  int nCallsPool;
  int nCallsHeur;
  int nCallsMWis1;
  int nCallsMWis2;
  int nCallsExact;
  double nTimePool;
  double nTimeHeur;
  double nTimeMwis1;
  double nTimeMwis2;
  double nTimeExact;

  Stats()
      : instance(""), solver(""), run(-1), nvertices(-1), nedges(-1), nA(-1),
        nB(-1), nvars(-1), ncons(-1), state(UNKNOWN), time(-1.0), nodes(0),
        nodesLeft(0), initSol(-1), initSolTime(-1.0), rootval(-1.0), lb(-1.0),
        ub(-1), gap(-1.0), ninfeas(0), ngcp(0), poolSize(-1), nColsPool(0),
        nColsHeur(0), nColsMwis1(0), nColsMwis2(0), nColsExact(0),
        nCallsPool(0), nCallsHeur(0), nCallsMWis1(0), nCallsMWis2(0),
        nCallsExact(0), nTimePool(0.0), nTimeHeur(0.0), nTimeMwis1(0.0),
        nTimeMwis2(0.0), nTimeExact(0.0) {}

  std::string get_state_as_str();
  void write_stats(std::ostream &file);
  void print_stats(std::ostream &file);
};

#endif // __STATS_HPP__