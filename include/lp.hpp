#ifndef _LP_HPP_
#define _LP_HPP_

#include "col.hpp"
#include "cplex_env.hpp"
#include "graph.hpp"
#include "params.hpp"
#include "pricing.hpp"
#include "stats.hpp"
extern "C" {
#include "color.h"
#include "color_private.h"
#include "mwis.h"
}

#include <chrono>
#include <iostream>

#define EPSILON 0.00001 // 10e-5

using Column = VertexVector;

class LP {

public:
  LP(Graph *graph, Params &params, Pool &pool, Graph &origGraph,
     std::ostream &log, bool isRoot = false);
  ~LP();

  // Optimize the linear relaxation by column generation
  [[nodiscard]] LP_STATE optimize(double timelimit, double ub, Stats &stats);

  // Get objective value (after calling optimize)
  [[nodiscard]] double get_obj_value() const { return objVal; };

  // Tell whether a feasible solution was found (after calling optimize)
  // Either found by the heuristic or the feasibility check
  [[nodiscard]] bool has_feas_sol() const {
    return initSol.get_n_colors() > 0;
  };

  // Get the cost of the feasible solution (after calling optimize)
  [[nodiscard]] size_t get_feas_value() const {
    return initSol.get_n_colors();
  };

  // Save the optimal solution (after calling optimize)
  void save_lp_solution(Stats &stats, Col &col);

  // Save the heuristic solution (after calling optimize)
  void save_heur_solution(Stats &stats, Col &col);

  // Branch (after calling optimize)
  void branch(std::vector<LP *> &branches);

private:
  // Input data
  GraphEnv in;
  // Parameters
  Params &params;
  // Vector of columns (stable sets)
  std::vector<Column> stables;
  // Vector of positive variables
  std::vector<int> posVars;
  // Branching variable
  Vertex branchVar;
  // Objective value
  double objVal;
  // Initial solution
  Col initSol;
  // Pool of columns
  Pool &pool;
  // Original graph
  Graph &origGraph;
  // Whether a dummy column was used to initialize the LP
  bool initializedWithDummy;
  // Log file
  std::ostream &log;
  // Is root the node?
  bool isRoot;
  // LP state
  LP_STATE state;

  // // Map from current vertices to original vertices
  // std::vector<Vertex> vertexMap;
  // // Map from original vertices to current vertices
  // std::vector<int> invVertexMap;

  // Heuristic solution of the DPCP instances at the current node
  void heuristic(HeurStats &stats, Params &params);

  // Check the feasibility of the DPCP instance at the current node
  void feasibility_check(Stats &stats, Params &params);

  // Initialize the linear relaxation with an initial set of columns
  void add_constraints_and_objective(CplexEnv &cenv);
  void add_initial_columns(CplexEnv &cenv);

  // Add a new column to the linear relaxation
  void add_column(CplexEnv &cenv, StableEnv &stab, bool addStable);

  // Check column validity
  bool check_column(StableEnv &stab);

  // Print column
  void print_column(StableEnv &stab);

  // Translate a stable set from the pool in terms of the vertices of the
  // current graph
  std::pair<bool, StableEnv> translate_stable_from_pool(StableEnv &stab,
                                                        IloNumArray &dualsA,
                                                        IloNumArray &dualsB);

  // Set CPLEX's parameters
  void set_parameters(CplexEnv &cenv, IloCplex &cplex);

  // Solve pricing problem
  int pricing(CplexEnv &cenv, PricingEnv &penv, Stats &stats,
              IloNumArray &dualsA, IloNumArray &dualsB, bool isRoot);
  int pricing_pool(CplexEnv &cenv, Stats &stats, IloNumArray &dualsA,
                   IloNumArray &dualsB, bool isRoot);
  int pricing_greedy(CplexEnv &cenv, PricingEnv &penv, Stats &stats,
                     IloNumArray &dualsA, IloNumArray &dualsB, bool isRoot);
  int pricing_P_Q_mwss(CplexEnv &cenv, PricingEnv &penv, Stats &stats,
                       IloNumArray &dualsA, IloNumArray &dualsB, bool isRoot);
  int pricing_P_mwss(CplexEnv &cenv, PricingEnv &penv, Stats &stats,
                     IloNumArray &dualsA, IloNumArray &dualsB, bool isRoot);
  int pricing_exact(CplexEnv &cenv, PricingEnv &penv, Stats &stats,
                    IloNumArray &dualsA, IloNumArray &dualsB, bool isRoot);

  // Exact solve of a GCP instance
  LP_STATE solve_GCP(Stats &stats, double timelimit, double ub);

  // Get branching variable
  Vertex get_branching_variable(const IloNumArray &values);
  Vertex get_branching_variable_FMS(const IloNumArray &values);
};

#endif // _LP_HPP_
