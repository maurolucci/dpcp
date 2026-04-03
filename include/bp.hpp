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
  explicit Node(LP&& lp, size_t depth = 0, size_t id = 0);

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

  void save(Col& sol);
  void save_heur(Col& sol);

  void branch(std::vector<Node>& sons);

 private:
  std::unique_ptr<LP> lp;
  size_t depth;
  size_t id;
};

class BP {
 public:
  BP(Params& params, std::ostream& log, Col& sol, double ub = DBL_MAX);

  Stats solve(Node root);

  // Methods for getting variables' value
  size_t get_nodes();
  double get_gap();
  double get_primal_bound();
  Stats& get_stats() { return stats; }

 private:
  std::list<Node> L;  // Priority queue

  Params& params;              // Input parameters
  Col& best_integer_solution;  // Current best integer solution
  double primal_bound;         // Primal bound
  size_t nodes;                // Number of processed nodes so far
  size_t nextNodeId;           // Next id to assign to a created tree node
  TimePoint start_t;           // B&P initial execution time
  TimePoint last_t;            // Used by log
  bool first_call;             // Used by log
  std::ostream& log;
  Stats stats;

  Stats return_stats(STATE state);
  LP_STATE push(Node node);
  Node& top();
  void pop();
  void update_primal_bound(double obj_value);
  double calculate_dual_bound();
  void show_stats();
};

#endif  // _BP_HPP_
