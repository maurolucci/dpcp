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

  // Heuristic options
  // heuristicRootNode: type of heuristic for the root node
  //    0: no heuristic
  //    1: greedy 1-step heuristic
  //    2: greedy 2-step heuristic
  //    3: semi-greedy 2-step heuristic (default)
  // heuristicOtherNodes: type of heuristic for other nodes
  //    0: no heuristic
  //    1: greedy 1-step heuristic
  //    2: greedy 2-step heuristic (default)
  //    3: semi-greedy 2-step heuristic
  // heuristic2stepVariant: variant of the 2-step heuristic
  //    2: DEG
  //    3: EDGE
  // heuristicSemigreedyAlpha: alpha parameter for semi-greedy heuristic
  // heuristicSemigreedyIter: number of iterations for semi-greedy heuristic
  int heuristicRootNode;
  int heuristicOtherNodes;
  size_t heuristic2stepVariant;
  double heuristicSemigreedyAlpha;
  size_t heuristicSemigreedyIter;

  // Feasibility check options
  // feasibilityRootNode: type of feasibility check for the root node
  //    0: no check
  //    1: enumerative
  //    2: ILP
  // feasibilityRootNodeTimeLimit: time limit for feasibility check (in
  // seconds). Only for ILP
  // feasibilityOtherNodes: type of feasibility check for other nodes
  //    0: no check
  //    1: enumerative
  //    2: ILP
  // feasibilityOtherNodesTimeLimit: time limit for feasibility check (in
  // seconds). Only for ILP
  int feasibilityRootNode;
  int feasibilityRootNodeTimeLimit;
  int feasibilityOtherNodes;
  int feasibilityOtherNodesTimeLimit;

  // Initialization options
  // inheritColumns: type of column inheritance from parent
  //     0: no inheritance
  //     1: inherit all columns
  //     2: inherit only positive columns
  // initializationBigWeight: weight of dummy column during initialization
  int inheritColumns;
  double initializationBigWeight;

  // Preprocessing options
  // preprocStep1: make each V_a a clique
  // preprocStep2: remove vertices for each a with |V_a| = 1
  // preprocStep3: vanish isolated vertices
  // preprocStep4: check if n = 1
  bool preprocStep1;
  bool preprocStep2;
  bool preprocStep3;
  bool preprocStep4;

  // Column generation options
  // usePool: use a pool of columns (not implemented yet)
  bool usePool;

  // Pricing options
  // pricingHeur1: use greedy heuristic for pricing
  // pricingHeur2: use P,Q-MWSSP heuristic for pricing
  // pricingHeur3: use P-MWSSP heuristic for pricing
  // pricingOrder: order of pricing
  //    1: pool -> greedy -> P,Q-MWSSP -> P-MWSSP
  //    2: pool -> greedy -> P-MWSSP -> P,Q-MWSSP
  // pricingHeur1MaxNCols: maximum number of columns to add with pricingHeur1
  // pricingExactTimeLimit: time limit for exact pricing (in seconds)
  bool pricingHeur1;
  bool pricingHeur2;
  bool pricingHeur3;
  int pricingOrder;
  size_t pricingHeur1MaxNCols;
  size_t pricingExactTimeLimit;
  // branchingFMS: use Furini-Malaguti-Santini branching rule
  bool branchingFMS;

  Params()
      : timeLimit(900), dfs(false), onlyRelaxation(false), heuristicRootNode(3),
        heuristicOtherNodes(2), heuristic2stepVariant(3),
        heuristicSemigreedyAlpha(0.2), heuristicSemigreedyIter(100),
        feasibilityRootNode(2), feasibilityRootNodeTimeLimit(300),
        feasibilityOtherNodes(2), feasibilityOtherNodesTimeLimit(60),
        inheritColumns(0), initializationBigWeight(1000.0), preprocStep1(true),
        preprocStep2(true), preprocStep3(true), preprocStep4(true),
        usePool(false), pricingHeur1(true), pricingHeur2(true),
        pricingHeur3(true), pricingOrder(1), pricingHeur1MaxNCols(1),
        pricingExactTimeLimit(300), branchingFMS(false){};

  std::string get_heur_name(int heur) {
    switch (heur) {
    case 0:
      return " (none)";
    case 1:
      return " (greedy 1-step)";
    case 2:
      return " (greedy 2-step)";
    case 3:
      return " (semi-greedy 2-step)";
    default:
      return " (unknown)";
    }
  };

  std::string get_heur_variant(int variant) {
    switch (variant) {
    case 0:
      return " (deg-real)";
    case 1:
      return " (deg-Q)";
    case 2:
      return " (deg-collapsed)";
    case 3:
      return " (edge)";
    default:
      return " (unknown)";
    }
  };

  std::string get_feas_name(int feas) {
    switch (feas) {
    case 0:
      return " (none)";
    case 1:
      return " (enumerative)";
    case 2:
      return " (ilp)";
    default:
      return " (unknown)";
    }
  };

  std::string get_inherit_name(int inherit) {
    switch (inherit) {
    case 0:
      return " (none)";
    case 1:
      return " (all)";
    case 2:
      return " (positive)";
    default:
      return " (unknown)";
    }
  };

  void print_params(std::ostream &out) {
    out << "*** Parameters ***:" << std::endl;
    out << "Time limit: " << timeLimit << " seconds" << std::endl;
    out << "DFS strategy: " << (dfs ? "enabled" : "disabled") << std::endl;
    out << "Only relaxation: " << (onlyRelaxation ? "enabled" : "disabled")
        << std::endl;
    out << "Heuristic root node: " << heuristicRootNode
        << get_heur_name(heuristicRootNode);
    if (heuristicRootNode == 2 || heuristicRootNode == 3)
      out << ", variant: " << heuristic2stepVariant
          << get_heur_variant(heuristic2stepVariant);
    if (heuristicRootNode == 3)
      out << ", alpha: " << heuristicSemigreedyAlpha
          << ", iterations: " << heuristicSemigreedyIter;
    out << std::endl;
    out << "Heuristic other nodes: " << heuristicOtherNodes
        << get_heur_name(heuristicOtherNodes);
    if (heuristicRootNode == 2 || heuristicRootNode == 3)
      out << ", variant: " << heuristic2stepVariant
          << get_heur_variant(heuristic2stepVariant);
    if (heuristicOtherNodes == 3)
      out << ", alpha: " << heuristicSemigreedyAlpha
          << ", iterations: " << heuristicSemigreedyIter;
    out << std::endl;
    out << "Feasibility root node: " << feasibilityRootNode
        << get_feas_name(feasibilityRootNode);
    if (feasibilityRootNode == 2)
      out << ", time limit: " << feasibilityRootNodeTimeLimit << " seconds";
    out << std::endl;
    out << "Feasibility other nodes: " << feasibilityOtherNodes
        << get_feas_name(feasibilityOtherNodes);
    if (feasibilityOtherNodes == 2)
      out << ", time limit: " << feasibilityOtherNodesTimeLimit << " seconds";
    out << std::endl;
    out << "Inherit columns: " << inheritColumns
        << get_inherit_name(inheritColumns) << std::endl;
    out << "Initialization big weight: " << initializationBigWeight
        << std::endl;
    out << "Preprocessing step 1 (clique V_a): "
        << (preprocStep1 ? "enabled" : "disabled") << std::endl;
    out << "Preprocessing step 2 (remove |V_a|=1): "
        << (preprocStep2 ? "enabled" : "disabled") << std::endl;
    out << "Preprocessing step 3 (vanish isolated): "
        << (preprocStep3 ? "enabled" : "disabled") << std::endl;
    out << "Preprocessing step 4 (check if n = 1): "
        << (preprocStep4 ? "enabled" : "disabled") << std::endl;
    out << "Use column pool: " << (usePool ? "enabled" : "disabled")
        << std::endl;
    out << "Pricing heuristic 1 (greedy): "
        << (pricingHeur1 ? "enabled" : "disabled")
        << ", max columns: " << pricingHeur1MaxNCols << std::endl;
    out << "Pricing heuristic 2 (P,Q-MWSSP): "
        << (pricingHeur2 ? "enabled" : "disabled") << std::endl;
    out << "Pricing heuristic 3 (P-MWSSP): "
        << (pricingHeur3 ? "enabled" : "disabled") << std::endl;
    out << "Pricing order: " << pricingOrder << std::endl;
    out << "Pricing exact time limit: " << pricingExactTimeLimit << " seconds"
        << std::endl;
  }
};

#endif // _PARAMS_HPP_