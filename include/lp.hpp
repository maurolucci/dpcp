#ifndef _LP_HPP_
#define _LP_HPP_

#include "col.hpp"
#include "cplex_env.hpp"
#include "graph.hpp"
#include "pricing.hpp"
#include "stats.hpp"
extern "C" {
#include "color.h"
#include "color_private.h"
#include "mwis.h"
}

#include <chrono>

#define EPSILON 0.00001     // 10e-5
#define MAXCOLSFROMPOOL 200 // Max number of columns added from pool per round
#define MAXCOLSFROMHEUR 200 // Max number of columns added from heur per round
#define MAXFAILSPOOL 100    // Max number of failed rounds for pool
#define MAXFAILSHEUR 100    // Max number of failed rounds for heur
#ifndef N_BRANCHES
#define N_BRANCHES 2
#endif

// typedef struct Column {
//   std::set<PSet> elements;
//   Column(int n_best, const nodepnt *best_sol);
// } Column;

class LP {

public:
  LP(const Graph &graph, Pool &pool, Graph &origGraph, Col *initSol = NULL);
  ~LP();

  // Optimize the linear relaxation by column generation
  [[nodiscard]] LP_STATE optimize(double timelimit);

  // Get objective value (after calling optimize)
  [[nodiscard]] double get_obj_value() const { return objVal; };

  // Get number of columns (after calling optimize)
  [[nodiscard]] size_t get_n_pool_columns() const { return nTotalPoolCols; };
  [[nodiscard]] size_t get_n_heur_columns() const { return nTotalHeurCols; };
  [[nodiscard]] size_t get_n_mwis2_columns() const { return nTotalMwis2Cols; };
  [[nodiscard]] size_t get_n_exact_columns() const { return nTotalExactCols; };

  // Save the optimal solution (after calling optimize)
  void save_solution(Col &col);

  // Branch (after calling optimize)
  void branch(std::vector<LP *> &branches);

private:
  // Input data
  GraphEnv in;
  // Vector of columns (stable sets)
  std::vector<std::vector<Vertex>> stables;
  // Vector of positive variables
  std::vector<int> posVars;
  // Branching variable
  Vertex branchVar;
  // Objective value
  double objVal;
  // Initial solution (only used in the root node)
  Col *initSol;
  // Pool of columns
  Pool &pool;
  // Original graph
  Graph &origGraph;
  // Map from current vertices to original vertices
  std::vector<Vertex> vertexMap;
  // Map from original vertices to current vertices
  std::vector<int> invVertexMap;
  size_t nTotalPoolCols;
  size_t nTotalHeurCols;
  size_t nTotalMwis2Cols;
  size_t nTotalExactCols;
  size_t nRemainingAttemptsPool;
  size_t nRemainingAttemptsHeur;
  size_t nRemainingAttemptsMwis2;

  // Initialize the linear relaxation with an initial set of columns
  void add_constraints(CplexEnv &cenv);
  void add_initial_columns(CplexEnv &cenv);

  // Add a new column to the linear relaxation
  void add_column(CplexEnv &cenv, StableEnv &stab);

  // Translate a stable set from the pool in terms of the vertices of the
  // current graph
  std::pair<bool, StableEnv> translate_stable_from_pool(StableEnv &stab,
                                                        IloNumArray &dualsA,
                                                        IloNumArray &dualsB);

  // Set CPLEX's parameters
  void set_parameters(CplexEnv &cenv, IloCplex &cplex);

  int pricing(CplexEnv &cenv, PricingEnv &penv, IloNumArray &dualsA,
              IloNumArray &dualsB);

  // // Compute vertex weights from dual values
  // // Return a boolean to indicate whether the weight of any vertex changed
  // sign
  // // from the previous iteration and the number of positive weights
  // std::pair<bool, size_t> get_weights(std::vector<double> &weights,
  //                                     IloNumArray &duals);

  // // Covert weights from double to int
  // int double2COLORNWT(COLORNWT nweights[], COLORNWT *scalef, size_t
  // nPosWeights,
  //                     const std::vector<double> &dbl_nweights);

  // Check if an stable from de pool is an entering column
  // Be careful, the stable is assume to be in terms of the original graph
  // and only coincides with the current graph at the root node
  bool check_stable(StableEnv &stab, IloNumArray &dualsA, IloNumArray &dualsB);

  // Exact solve of a GCP instance
  LP_STATE solve_GCP(double timelimit);

  // Get branching variable
  size_t get_branching_variable(const IloNumArray &values);
};

#endif // _LP_HPP_
