#ifndef _PARAMS_HPP_
#define _PARAMS_HPP_

#include <string>

struct Params {
  // General options
  // timeLimit: overall time limit for the algorithm (in seconds)
  // dfs: use depth-first strategy in the B&P tree
  // onlyRelaxation: solve only the root node
  size_t timeLimit;
  bool dfs;
  bool onlyRelaxation;
  // Initialization options
  // initializationBigWeight: weight of dummy column during initialization
  // initializationUseHeur: use heuristic to find an initial solution
  double initializationBigWeight;
  bool initializationUseHeur;
  // Preprocessing options
  // preprocStep1: make each V_a a clique
  // preprocStep2: remove vertices for each a with |V_a| = 1
  // preprocStep3: vanish isolated vertices
  bool preprocStep1;
  bool preprocStep2;
  bool preprocStep3;
  // Pool of columns
  bool usePool;
  // Pricing options
  // pricingHeur1: use greedy heuristic for pricing
  // pricingHeur2: use MWSSP I heuristic for pricing
  // pricingHeur3: use MWSSP II heuristic for pricing
  // pricingOrder: order of pricing
  //    1: pool -> greedy -> MWSSP I -> MWSSP II
  //    2: pool -> greedy -> MWSSP II -> MWSSP I
  // pricingHeur1MaxNCols: maximum number of columns to add with pricingHeur1
  // pricingExactTimeLimit: time limit for exact pricing (in seconds)
  bool pricingHeur1;
  bool pricingHeur2;
  bool pricingHeur3;
  int pricingOrder;
  size_t pricingHeur1MaxNCols;
  size_t pricingExactTimeLimit;
  Params()
      : timeLimit(900), dfs(false), onlyRelaxation(false),
        initializationBigWeight(1000.0), initializationUseHeur(false),
        preprocStep1(true), preprocStep2(true), preprocStep3(true),
        usePool(false), pricingHeur1(true), pricingHeur2(true),
        pricingHeur3(true), pricingOrder(1), pricingHeur1MaxNCols(1),
        pricingExactTimeLimit(300){};
};

#endif // _PARAMS_HPP_