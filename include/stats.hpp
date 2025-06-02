#ifndef __STATS_HPP__
#define __STATS_HPP__

#include <iostream>
#include <limits>
#include <string>

enum PRICING_STATE {
  PRICING_SOLUTION,
  PRICING_NO_SOLUTION,
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
  UNKNOWN,
};

struct Stats {
  int nvars;
  int ncons;
  STATE state;
  double time;
  int nodes;
  int initSol;
  double lb;
  int ub;
  double gap;
  int poolSize;
  int ncolsPool;
  int ncolsHeur;
  int ncolsMwis2;
  int ncolsExact;

  Stats()
      : nvars(-1), ncons(-1), state(UNKNOWN), time(-1.0), nodes(-1),
        initSol(-1), lb(-1.0), ub(-1), gap(-1.0), poolSize(-1), ncolsPool(-1),
        ncolsHeur(-1), ncolsMwis2(-1), ncolsExact(-1){};

  std::string get_state_as_str();
  void write_stats(std::ostream &file);
  void print_stats(std::ostream &file);
};

#endif // __STATS_HPP__