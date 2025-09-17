#ifndef __PRICING_HPP__
#define __PRICING_HPP__

#include "graph.hpp"
#include "stats.hpp"
extern "C" {
#include "color.h"
#include "color_private.h"
#include "mwis.h"
}

#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>
#include <random>

#define THRESHOLD 1.1           // Threshold for early stop
#define PRICING_EPSILON 0.00001 // 10e-5

// This is the class implementing the generic callback interface.
class ThresholdCallback : public IloCplex::Callback::Function {

private:
  GraphEnv &in;
  StableEnv &stab;

  // Variables
  IloNumVarArray &y;
  IloNumVarArray &w;

public:
  // Constructor with data.
  ThresholdCallback(GraphEnv &in, StableEnv &stab, IloNumVarArray &y,
                    IloNumVarArray &w)
      : in(in), stab(stab), y(y), w(w){};

  inline void check_thresolhd(const IloCplex::Callback::Context &context);

  // This is the function that we have to implement and that CPLEX will call
  // during the solution process at the places that we asked for.
  virtual void invoke(const IloCplex::Callback::Context &context) ILO_OVERRIDE;

  /// Destructor
  ~ThresholdCallback(){};
};

class PricingEnv {

private:
  GraphEnv &in;
  StableEnv stab;
  double exactTimeLimit;

  // CPLEX variables
  IloEnv cxenv;
  IloModel cxmodel;
  IloObjective cxobj;
  IloNumVarArray y;
  IloNumVarArray w;
  IloConstraintArray cxcons;
  IloCplex cplex;

  // Callback variables
  ThresholdCallback cb;
  CPXLONG contextMask;

  // MWIS variables
  MWISenv *mwis_env;
  COLORNWT *mwis_pi;
  int ecount;
  int *elist;

  void exact_init();
  void mwis_init();

  // std::pair<bool, size_t> get_weights(std::vector<double> &weights,
  //                                    IloNumArray &dualsA, IloNumArray
  //                                    &dualsB);

  int double2COLORNWT(COLORNWT nweights[], COLORNWT *scalef,
                      const std::vector<double> &dbl_weights);

public:
  PricingEnv(GraphEnv &in, double exactTimeLimit);
  ~PricingEnv();

  std::pair<StableEnv, PRICING_STATE> heur_solve(IloNumArray &dualsA,
                                                 IloNumArray &dualsB);

  std::list<std::pair<StableEnv, PRICING_STATE>>
  mwis1_solve(IloNumArray &dualsA, IloNumArray &dualsB);
  std::pair<StableEnv, PRICING_STATE> mwis2_solve(IloNumArray &dualsA,
                                                  IloNumArray &dualsB);

  std::pair<StableEnv, PRICING_STATE> exact_solve(IloNumArray &dualsA,
                                                  IloNumArray &dualsB);
};

#endif // __PRICING_HPP__