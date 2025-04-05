#ifndef __PRICING_HPP__
#define __PRICING_HPP__

#include "graph.hpp"
#include "stats.hpp"

#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>
#include <random>

#define TIMELIMIT 300.0         // = 5 minutes
#define THRESHOLD 1.1           // Threshold for early stop
#define PRICING_EPSILON 0.00001 // 10e-5

typedef std::mt19937 rng_type;

// This is the class implementing the generic callback interface.
class ThresholdCallback : public IloCplex::Callback::Function {

private:
  GraphEnv &in;
  StableEnv &stab;

  // Variables
  IloNumVarArray &y;
  IloNumVarArray &wa;
  IloNumVarArray &wb;

public:
  // Constructor with data.
  ThresholdCallback(GraphEnv &in, StableEnv &stab, IloNumVarArray &y,
                    IloNumVarArray &wa, IloNumVarArray &wb)
      : in(in), stab(stab), y(y), wa(wa), wb(wb) {};

  inline void check_thresolhd(const IloCplex::Callback::Context &context);

  // This is the function that we have to implement and that CPLEX will call
  // during the solution process at the places that we asked for.
  virtual void invoke(const IloCplex::Callback::Context &context) ILO_OVERRIDE;

  /// Destructor
  ~ThresholdCallback() {};
};

class PricingEnv {

private:
  GraphEnv &in;
  StableEnv stab;

  // CPLEX variables
  IloEnv cxenv;
  IloModel cxmodel;
  IloObjective cxobj;
  IloNumVarArray y;
  IloNumVarArray wa;
  IloNumVarArray wb;
  IloConstraintArray cxcons;
  IloCplex cplex;

  // Callback variables
  ThresholdCallback cb;
  CPXLONG contextMask;

  // Random
  rng_type rng;

  // Heuristic pricing variables
  std::list<std::tuple<double, size_t, TypeA, TypeB>> heurCandidates;

  void exact_init();

  void try_stable_add(double weight, size_t v, TypeA a, TypeB b,
                      IloNumArray &dualsA, std::vector<bool> &used,
                      std::set<TypeB> &bs);

public:
  PricingEnv(GraphEnv &in);
  ~PricingEnv();

  void heur_init(IloNumArray &dualsA, IloNumArray &dualsB);

  std::pair<StableEnv, PRICING_STATE>
  heur_solve(IloNumArray &dualsA, IloNumArray &dualsB, TypeA first);

  std::pair<StableEnv, PRICING_STATE> exact_solve(IloNumArray &dualsA,
                                                  IloNumArray &dualsB);
};

#endif // __PRICING_HPP__