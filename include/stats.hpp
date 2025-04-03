#ifndef __STATS_HPP__
#define __STATS_HPP__

#include <iostream>
#include <string>

enum PRICING_STATE {
  PRICING_OPTIMAL,
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
  LP_MEM_EXCEEDED,
};

enum STATE {
  OPTIMAL,
  FEASIBLE,
  INFEASIBLE,
  TIME_EXCEEDED,
  MEM_EXCEEDED,
  NODE_TIME_EXCEEDED,
  NODE_MEM_EXCEEDED,
  UNKNOWN,
};

struct Stats {
  size_t nvars;
  size_t ncons;
  STATE state;
  double time;
  size_t nodes;
  double lb;
  double ub;
  double gap;

  std::string get_state_as_str();
  void write_stats(std::ostream &file);
  void print_stats(std::ostream &file);
};

#endif // __STATS_HPP__