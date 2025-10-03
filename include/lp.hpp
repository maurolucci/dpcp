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
  LP(const Graph &graph, Params &params, Pool &pool, Graph &origGraph,
     std::ostream &log, bool isRoot = false);
  ~LP();

  // Optimize the linear relaxation by column generation
  [[nodiscard]] LP_STATE optimize(double timelimit, Stats &stats);

  // Get objective value (after calling optimize)
  [[nodiscard]] double get_obj_value() const { return objVal; };

  // Apply heuristic
  [[nodiscard]] bool solve_heur(double &obj_value, double &time, bool isRoot);

  // Save the optimal solution (after calling optimize)
  void save_lp_solution(Col &col);

  // Save the heuristic solution (after calling solve_heur)
  void save_heur_solution(Col &col);

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

  // // Map from current vertices to original vertices
  // std::vector<Vertex> vertexMap;
  // // Map from original vertices to current vertices
  // std::vector<int> invVertexMap;

  void solve_heur_aux(Stats &stats, int heur, int iters);

  // Initialize the linear relaxation with an initial set of columns
  void add_constraints_and_objective(CplexEnv &cenv);
  void add_initial_columns(CplexEnv &cenv);

  // Add a new column to the linear relaxation
  void add_column(CplexEnv &cenv, StableEnv &stab);

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
              IloNumArray &dualsA, IloNumArray &dualsB);
  int pricing_pool(CplexEnv &cenv, Stats &stats, IloNumArray &dualsA,
                   IloNumArray &dualsB);
  int pricing_greedy(CplexEnv &cenv, PricingEnv &penv, Stats &stats,
                     IloNumArray &dualsA, IloNumArray &dualsB);
  int pricing_mwss1(CplexEnv &cenv, PricingEnv &penv, Stats &stats,
                    IloNumArray &dualsA, IloNumArray &dualsB);
  int pricing_mwss2(CplexEnv &cenv, PricingEnv &penv, Stats &stats,
                    IloNumArray &dualsA, IloNumArray &dualsB);
  int pricing_exact(CplexEnv &cenv, PricingEnv &penv, Stats &stats,
                    IloNumArray &dualsA, IloNumArray &dualsB);

  // Exact solve of a GCP instance
  LP_STATE solve_GCP(double timelimit);

  // Get branching variable
  Vertex get_branching_variable(const IloNumArray &values);
};

#endif // _LP_HPP_
