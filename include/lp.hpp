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

#define EPSILON 0.00001  // 10e-5

using Column = StableEnv;
using Pool = std::vector<Column>;

class LP {
 public:
  LP(DPCPInst dpcp, Pool pool, const DPCPInst& origDpcp, Params& params,
     Stats& stats, std::ostream& log, bool isRoot = false);
  ~LP();

  // Optimize the linear relaxation by column generation
  [[nodiscard]] LP_STATE solve(double timelimit, double ub);

  // Getters
  [[nodiscard]] double get_lower_bound() const { return objVal; };
  [[nodiscard]] size_t get_upper_bound() const {
    return coloring.get_n_colors();
  };
  [[nodiscard]] bool has_heur_solution() const {
    return coloring.get_n_colors() > 0;
  };
  Col get_heur_solution();
  Col get_lp_solution();

  // Branch
  std::vector<LP> branch();

 private:
  DPCPInst dpcp;  // DPCP instance at the current node
  Pool pool;      // Pool of stable sets at the current node (herited from the
                  // parent)
  const DPCPInst& origDpcp;  // Reference to the original DPCP instance
  Params& params;      // Reference to the parameters of the algorithm
  Stats& stats;  // Reference to the stats object to update during the algorithm
  std::ostream& log;  // Reference to the log stream

  bool isRoot;                // Is the current node the root node?
  double objVal;              // Objective value of the current LP
  LP_STATE state;             // State of the current LP
  bool initializedWithDummy;  // Was the current LP initialized with a dummy
                              // column?

  std::vector<Column>
      stables;  // Stable sets corresponding to the columns in the current LP
  std::vector<size_t> posVars;  // Indices of columns whose variables have
                                // positive value in the current LP solution
  Vertex branchingVertex;  // Vertex selected for branching at the current node

  Col coloring;  // Feasible coloring for the current node

  // Exact solve of a GCP instance
  LP_STATE gcp_solve(double timelimit, double ub);

  // Heuristic solution of the DPCP instances at the current node
  void heuristic_solve();

  // Feasibility check of the DPCP instance at the current node
  bool feasibility_solve();

  // Set CPLEX's parameters
  void set_parameters(CplexEnv& cenv, IloCplex& cplex);

  // Initialize LP
  void add_constraints_and_objective(CplexEnv& cenv);
  void add_initial_columns(CplexEnv& cenv);

  // Pricing methods
  int pricing(CplexEnv& cenv, PricingEnv& penv, IloNumArray& dualsP,
              IloNumArray& dualsQ);
  int pricing_pool(CplexEnv& cenv, IloNumArray& dualsP, IloNumArray& dualsQ);
  int pricing_greedy(CplexEnv& cenv, PricingEnv& penv, IloNumArray& dualsP,
                     IloNumArray& dualsQ, bool enabled);
  int pricing_P_Q_mwss(CplexEnv& cenv, PricingEnv& penv, IloNumArray& dualsP,
                       IloNumArray& dualsQ, bool enabled);
  int pricing_P_mwss(CplexEnv& cenv, PricingEnv& penv, IloNumArray& dualsP,
                     IloNumArray& dualsQ, bool enabled);
  int pricing_exact(CplexEnv& cenv, PricingEnv& penv, IloNumArray& dualsP,
                    IloNumArray& dualsQ);

  // Add a new column to the linear relaxation
  void add_column(CplexEnv& cenv, Column& stab);

  // Check column validity
  bool check_column(Column& stab);

  // Print column
  void print_column(Column& stab);

  // Get branching variable
  Vertex get_branching_variable_LNTT(const IloNumArray& values);
  Vertex get_branching_variable_FMS(const IloNumArray& values);
};

#endif  // _LP_HPP_
