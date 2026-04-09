#ifndef _BP_HPP_
#define _BP_HPP_

#include <cfloat>
#include <chrono>
#include <iostream>
#include <list>
#include <memory>
#include <vector>

#include "col.hpp"
#include "lp.hpp"
#include "stats.hpp"

#define EPSILON_BP 0.001  // For doing ceil(x - EPSILON_BP) during prunning

using ClockType = std::chrono::high_resolution_clock;
using TimePoint = ClockType::time_point;

class Node {
 public:
  Node(const DPCPInst& origDpcp, Params& params, Stats& stats,
  std::ostream& log, std::ostream& debugLog, std::ostream& colLog,
  bool isRoot = false);
  explicit Node(Node& parent, BRANCH_NODE branchNode, size_t depth = 0,
                size_t id = 0);

  Node(const Node&) = delete;
  Node& operator=(const Node&) = delete;
  Node(Node&&) noexcept = default;
  Node& operator=(Node&&) noexcept = default;

  double get_obj_value() const;
  size_t get_depth() const { return depth; }
  size_t get_id() const { return id; }
  void set_id(size_t nodeId) { id = nodeId; }
  bool operator>(const Node& n) const;

  LP_STATE solve(double timelimit, double ub);

  bool feas_sol() const;
  size_t feas_value() const;
  LP_INTEGER_SOURCE integer_source() const;

  void save(Col& sol);
  void save_heur(Col& sol);
  LP& get_lp() { return lp; }

  void branch(std::vector<std::unique_ptr<Node>>& sons);

 private:
  LP lp;
  size_t depth;
  size_t id;
};

class BP {
 public:
  BP(Params& params, std::ostream& log, std::ostream& debugLog,
    std::ostream& colLog, Col& sol, double ub = DBL_MAX);

  Stats solve(DPCPInst& dpcp);

  // Methods for getting variables' value
  size_t get_nodes();
  double get_gap();
  double get_primal_bound();
  Stats& get_stats() { return stats; }

 private:
  std::list<std::unique_ptr<Node>> L;  // Priority queue

  Params& params;              // Input parameters
  Col& best_integer_solution;  // Current best integer solution
  double primal_bound;         // Primal bound
  size_t nodes;                // Number of processed nodes so far
  size_t nextNodeId;           // Next id to assign to a created tree node
  TimePoint start_t;           // B&P initial execution time
  TimePoint last_t;            // Used by log
  bool first_call;             // Used by log
  std::ostream& log;
  std::ostream& debugLog;
  std::ostream& colLog;
  Stats stats;

  Stats return_stats(STATE state);
  LP_STATE push(std::unique_ptr<Node> node);
  Node& top();
  void pop();
  void update_primal_bound(double obj_value);
  double calculate_dual_bound();
  void show_stats();
};

#endif  // _BP_HPP_
