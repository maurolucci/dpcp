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

using ClockType = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::_V2::system_clock::time_point;

class Node {

public:
  Node(LP *lp) : lp(lp) {}

  ~Node() { delete lp; }

  double get_obj_value() const { return lp->get_obj_value(); }

  bool operator>(const Node &n) const {
    return (get_obj_value() > n.get_obj_value());
  }

  LP_STATE solve(double timelimit, double ub, Stats &stats) {
    LP_STATE state = lp->optimize(timelimit, ub, stats);
    return state;
  }

  bool feas_sol() const { return lp->has_feas_sol(); }
  size_t feas_value() const { return lp->get_feas_value(); }

  template <class Solution> void save(Solution &sol) {
    lp->save_lp_solution(sol);
  }

  template <class Solution> void save_heur(Solution &sol) {
    lp->save_heur_solution(sol);
  }

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
  BP(Params &params, std::ostream &log, Solution &sol, double ub = DBL_MAX)
      : params(params), best_integer_solution(sol), primal_bound(ub), nodes(0),
        log(log), stats() {}

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
    case LP_INIT_FAIL:
      return return_stats(INIT_FAIL);
    default:
      break;
    }

    if (!params.onlyRelaxation) {
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
          case LP_INIT_FAIL:
            delete node;
            return return_stats(INIT_FAIL);
          default:
            break;
          }
        }

        delete node;
      }
    }

    if (primal_bound == DBL_MAX) // Infeasibility case:
      return return_stats(INFEASIBLE);

    if (params.onlyRelaxation && !L.empty())
      return return_stats(FEASIBLE);

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

  Params &params;                  // Input parameters
  Solution &best_integer_solution; // Current best integer solution
  double primal_bound; // Primal bound (given by the best integer solution)
  size_t
      nodes; // Number of processed nodes so far. A node is considered processed
  // if its relaxation has been solved
  TimePoint start_t; // B&P initial execution time
  TimePoint last_t;  // used by log
  bool first_call;   // used by log
  int opt_flag;      // Optimality flag
  std::ostream &log;
  Stats stats;

  Stats return_stats(STATE state) {
    stats.state = state;
    stats.time =
        std::chrono::duration<double>(ClockType::now() - start_t).count();
    if (stats.time > params.timeLimit) {
      stats.time = params.timeLimit;
      stats.state = primal_bound == DBL_MAX ? TIME_EXCEEDED : FEASIBLE;
    }
    stats.nodes = static_cast<int>(nodes);
    stats.nodesLeft = static_cast<int>(L.size());
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

    return stats;
  }

  LP_STATE push(Node *node, bool root = false) {

    nodes++;
    double obj_value;

    // Solve the linear relaxation of the node and prune if possible
    double elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                         ClockType::now() - start_t)
                         .count();
    LP_STATE state =
        node->solve(params.timeLimit - elapsed, primal_bound, stats);

    // If a better feasible solution was found, save it
    if (node->feas_sol()) {
      obj_value = node->feas_value();
      if (obj_value < primal_bound) {
        node->save_heur(best_integer_solution);
        log << "New best integer solution found by heuristic with value: "
            << obj_value << std::endl;
        update_primal_bound(obj_value);
      }
    }

    switch (state) {

    case LP_INTEGER:
      obj_value = node->get_obj_value();
      if (obj_value < primal_bound) {
        node->save(best_integer_solution);
        log << "New best integer solution found by LR with value: " << obj_value
            << std::endl;
        update_primal_bound(obj_value);
      }
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

    if (params.dfs) { // DFS strategy
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
      if (*node > **it) {
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

  void update_primal_bound(double obj_value) {

    // Update best value
    primal_bound = obj_value;

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
      if (params.dfs) {
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
