#ifndef _BP_HPP_
#define _BP_HPP_

#include "bp.hpp"
#include "lp.hpp"
#include "stats.hpp"

#include <cfloat>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iterator>

#define EPSILON_BP 0.001 // For doing ceil(x - EPSILON_BP) during prunning

// #define ONLY_RELAXATION
#ifndef MAXTIME
#define MAXTIME 900 // In seconds
#endif

using ClockType = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::_V2::system_clock::time_point;

class Node {

public:
  Node(LP *lp) : lp(lp) {}

  ~Node() { delete lp; }

  double get_obj_value() const { return lp->get_obj_value(); }

  bool operator<(const Node &n) const {
    return (get_obj_value() < n.get_obj_value());
  }

  LP_STATE solve(double timelimit, Stats &stats) {
    LP_STATE state = lp->optimize(timelimit, stats);
    return state;
  }

  template <class Solution> void save(Solution &sol) { lp->save_solution(sol); }

  void branch(std::vector<Node *> &sons) {

    std::vector<LP *> lps;
    lp->branch(lps);

    sons.reserve(lps.size());
    for (auto x : lps)
      sons.push_back(new Node(x));

    return;
  }

private:
  LP *lp;
};

template <class Solution> class BP {

public:
  BP(Solution &sol, std::ostream &log, bool DFS = false)
      : best_integer_solution(sol), primal_bound(DBL_MAX), nodes(0), DFS(DFS),
        log(log), stats(), initSolValue(-1.0) {}

  void set_initial_solution(Solution &initSol, double initValue) {
    initSolValue = primal_bound = initValue;
    best_integer_solution = initSol;
  }

  Stats solve(Node *root) {

    start_t = ClockType::now();
    last_t = start_t;
    first_call = true;

    // Push root note (and solve initial LR)
    switch (push(root, true)) {
    case LP_TIME_EXCEEDED:
      return return_stats(TIME_EXCEEDED_LP);
    case LP_TIME_EXCEEDED_PR:
      return return_stats(TIME_EXCEEDED_PR);
    case LP_MEM_EXCEEDED:
      return return_stats(MEM_EXCEEDED_LP);
    case LP_MEM_EXCEEDED_PR:
      return return_stats(MEM_EXCEEDED_PR);
    default:
      break;
    }

#ifndef ONLY_RELAXATION
    while (!L.empty()) {

      // Pop next node
      Node *node = top();
      show_stats(*node); // First show_stats, then pop
      pop();

      // Re-try to prune by bound, since primal_bound could have been improved
      if (ceil(node->get_obj_value() - EPSILON_BP) >= primal_bound) {
        delete node;
        continue;
      }

      // Branch
      std::vector<Node *> sons;
      node->branch(sons);

      // Push sons (and solve initial LR)
      for (auto n : sons) {

        switch (push(n)) {
        case LP_TIME_EXCEEDED:
          delete node;
          return return_stats(TIME_EXCEEDED_LP);
        case LP_TIME_EXCEEDED_PR:
          delete node;
          return return_stats(TIME_EXCEEDED_PR);
        case LP_MEM_EXCEEDED:
          delete node;
          return return_stats(MEM_EXCEEDED_LP);
        case LP_MEM_EXCEEDED_PR:
          delete node;
          return return_stats(MEM_EXCEEDED_PR);
        default:
          break;
        }
      }

      delete node;
    }
#endif

    if (primal_bound == DBL_MAX) // Infeasibility case:
      return return_stats(INFEASIBLE);

#ifdef ONLY_RELAXATION
    if (!L.empty())
      return return_stats(FEASIBLE);
#endif

    // Optimality case:
    return return_stats(OPTIMAL);
  }

  // Methods for getting variables' value
  size_t get_nodes() { return nodes; }

  double get_gap() {
    double _dual_bound = calculate_dual_bound();
    double gap = DBL_MAX;
    if (_dual_bound != -DBL_MAX)
      gap = abs(_dual_bound - primal_bound) /
            (0.0000000001 + abs(primal_bound)) * 100;
    return gap;
  }

  double get_primal_bound() { return primal_bound; }

  int get_opt_flag() { return opt_flag; }

private:
  std::list<Node *> L; // Priority queue

  Solution &best_integer_solution; // Current best integer solution
  double primal_bound; // Primal bound (given by the best integer solution)
  size_t
      nodes; // Number of processed nodes so far. A node is considered processed
  // if its relaxation has been solved
  TimePoint start_t; // B&P initial execution time
  TimePoint last_t;  // used by log
  bool first_call;   // used by log
  int opt_flag;      // Optimality flag
  bool DFS;
  std::ostream &log;
  Stats stats;
  double initSolValue; // Value of the provided initial solution
  double rootval;

  Stats return_stats(STATE state) {
    stats.state = state;
    stats.time =
        std::chrono::duration<double>(ClockType::now() - start_t).count();
    if (stats.time > MAXTIME) {
      stats.time = MAXTIME;
      stats.state = primal_bound == DBL_MAX ? TIME_EXCEEDED : FEASIBLE;
    }
    stats.nodes = static_cast<int>(nodes);
    stats.rootval = rootval;
    stats.ub = static_cast<int>(primal_bound + 0.5);
    if (state == OPTIMAL) {
      stats.lb = primal_bound;
      stats.gap = 0.0;
    } else if (state == INFEASIBLE) {
      ;
    } else {
      stats.lb = calculate_dual_bound();
      stats.gap = get_gap() / 100;
    }
    stats.initSol = static_cast<int>(initSolValue + 0.5);

    return stats;
  }

  LP_STATE push(Node *node, bool root = false) {

    nodes++;
    double obj_value;

    // Solve the linear relaxation of the node and prune if possible
    double elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                         ClockType::now() - start_t)
                         .count();
    LP_STATE state = node->solve(MAXTIME - elapsed, stats);

    if (root)
      rootval = node->get_obj_value();

    switch (state) {

    case LP_INTEGER:
      obj_value = node->get_obj_value();
      if (obj_value < primal_bound)
        update_primal_bound(*node);
      delete node; // Prune by optimality
      return state;

    case LP_FRACTIONAL:
      obj_value = node->get_obj_value();
      if (ceil(obj_value - EPSILON_BP) >= primal_bound) {
        delete node; // Prune by optimality
        return state;
      }
      break; // Do not prune

    default:
      delete node; // Prune by infeasibility or mem/time limit
      return state;
    }

    // Continuation for fractional nodes
    // Place the node in the list according to its priority

    if (DFS) { // DFS strategy
      L.push_back(node);
      return LP_FRACTIONAL;
    }

    // Otherwise, best-bound strategy

    // list is empty
    if (L.empty()) {
      L.push_back(node);
      return state;
    }

    // list is not empty
    for (auto it = L.begin(); it != L.end(); ++it)
      if (*node < **it) {
        L.insert(it, node);
        return state;
      }

    // Otherwise, push back
    L.push_back(node);
    return state;
  }

  Node *top() { return L.back(); }

  void pop() {
    L.pop_back();
    return;
  }

  void update_primal_bound(Node &node) {

    // Update best integer solution and value
    node.save(best_integer_solution);
    primal_bound = node.get_obj_value();

    // Prune if possible
    for (auto it = L.begin(); it != L.end();) {
      if ((*it)->get_obj_value() >= primal_bound) {
        delete *it;
        it = L.erase(it);
      } else
        ++it;
    }
  }

  double calculate_dual_bound() {

    double _dual_bound = DBL_MAX; // minimum objective value of unpruned nodes
    if (!L.empty()) {
      if (DFS) {
        // Traverse the list and search for the minimum objective value
        for (auto it = L.begin(); it != L.end(); ++it)
          if ((*it)->get_obj_value() < _dual_bound)
            _dual_bound = (*it)->get_obj_value();
      } else {
        _dual_bound = L.front()->get_obj_value();
      }
    }

    if (_dual_bound == DBL_MAX) {
      return -DBL_MAX;
    } else {
      return _dual_bound;
    }
  }

  void show_stats(Node &node) {

    auto now_t = ClockType::now();

    if (first_call)
      first_call = false;
    else {
      if (std::chrono::duration<double>(now_t - last_t).count() < 10.0)
        return;
      last_t = now_t;
    }

    double _dual_bound =
        calculate_dual_bound(); // (note: it is time consuming when DFS is used)

    log << std::fixed << std::setprecision(3);

    log << "  Best obj value  = " << _dual_bound << "\t Best int = ";
    if (primal_bound == DBL_MAX)
      log << "inf\t Gap = ---";
    else
      log << (int)(EPSILON + primal_bound) << "\t Gap = "
          << (primal_bound - _dual_bound) / (EPSILON + primal_bound) * 100
          << "%";
    log << "\t Nodes: processed = " << nodes << ", left = " << L.size()
        << "\t time = "
        << std::chrono::duration<double>(now_t - start_t).count() << std::endl;
  }
};

#endif // _BP_HPP_
