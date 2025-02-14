#ifndef _LP_HPP_
#define _LP_HPP_

#include "col.hpp"
#include "graph.hpp"
extern "C" {
#include "mwis.h"
}

#include <chrono>
#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>
#include <set>
#include <vector>

#define EPSILON 0.00001
#define MAXTIME 300.0
#define THRESHOLD 0.1
// #define IL_STD 1 // CPLEX lo pide
#ifndef N_BRANCHES
#define N_BRANCHES 2
#endif

enum LP_STATE { INFEASIBLE, INTEGER, FRACTIONAL, TIME_OR_MEM_LIMIT };

// typedef struct Column {
//   std::set<PSet> elements;
//   Column(int n_best, const nodepnt *best_sol);
// } Column;

class LP {

public:
  // Constructors
  LP(const HyperGraph &hg, const ConflictGraph &cg, double start_t);
  LP(const HyperGraph &hg, const ConflictGraph &cg, int vcount, CGVertex *vlist,
     int ecount, int *elist, double start_t);

  // Destructor
  ~LP();

  // Optimize LP
  LP_STATE optimize();

  // Save optimal solution
  void save_solution(Coloring &coloring);

  // Get objective value
  [[nodiscard]] inline double get_obj_value() const { return obj_value; };

  // Get number of columns
  [[nodiscard]] inline int get_n_columns() const { return Xvars.getSize(); };

  // Branch
  void branch(std::vector<LP *> &branches, CGVertex cgv);

private:
  double start_t;       // Starting time of B&P
  IloEnv Xenv;          // CPLEX environment structure
  IloModel Xmodel;      // CPLEX model
  IloObjective Xobj;    // CPLEX objective function
  IloNumVarArray Xvars; // CPLEX variables
  IloRangeArray Xrestr; // CPLEX constraints
  // IloArray<Column> vars; // Internal representation of CPLEX's columns
  IloNumArray values; // Value of the optimal solution
  IloNum obj_value;   // Objective value of the optimal solution

  std::list<int> pos_vars; // List with the indexes of the positive variables
  int most_fract_var;      // Index of the most fractional var with at least two
                           // vertices in the stable set

  const HyperGraph &hg;    // Original hypergraph
  const ConflictGraph &cg; // Original conflict graph
  size_t nrows;            // Number of constraints (= number of dual variables)
  int vcount;
  CGVertex *vlist;
  int ecount;
  int *elist;

  // Initialize LP
  void initialize();

  int double2COLORNWT(COLORNWT nweights[], COLORNWT *scalef,
                      const IloNumArray dbl_nweights);

  // Fill the LP with a feasible set of columns
  void fill_initial_columns();

  void add_column(int count, const int *members);

  void set_parameters(IloCplex &cplex);
};

#endif // _LP_HPP_
