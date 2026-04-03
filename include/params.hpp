#ifndef _PARAMS_HPP_
#define _PARAMS_HPP_

#include <string>

struct Params {
  // General options
  // timeLimit: overall time limit for the algorithm (in seconds)
  // treeSearch: tree search strategy in B&P
  //   1: best bound
  //   2: depth-first search
  // onlyRelaxation: solve only the root node
  // verbose: verbosity level (0: quiet, 1: low, 2: detailed)
  size_t timeLimit;
  int treeSearch;
  bool onlyRelaxation;
  int verbose;

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
  // preprocessing: run all preprocessing steps (clique V_a, remove |V_a|=1,
  //   vanish isolated vertices, check if n=1)
  bool preprocessing;

  // Pricing options
  // pricingMethod: pricing strategy selector
  //    0: ILP
  //    1: greedy -> ILP
  //    2: greedy -> P,Q-MWSSP -> ILP
  //    3: greedy -> P-MWSSP -> ILP
  //    4: greedy -> P,Q-MWSSP -> P-MWSSP -> ILP
  //    5: greedy -> P-MWSSP -> P,Q-MWSSP -> ILP
  //    6: automatic (depends on graph density)
  // pricingHeur1MaxNCols: maximum number of columns to add with greedy
  // pricing
  // pricingMaxColsPerIter: maximum number of columns added per pricing call
  //   (applies to pool and greedy; best k selected via partial sort)
  // pricingExactTimeLimit: time limit for exact pricing (in seconds)
  int pricingMethod;
  double pricingHeur1Alpha;
  size_t pricingHeur1MaxNCols;
  size_t pricingMaxColsPerIter;
  size_t pricingExactTimeLimit;
  // branchingVariable: branching variable rule
  //   1: Furini-Malaguti-Santini (FMS)
  //   2: Lucci-Nasini-Tolomei-Torres (LNTT)
  int branchingVariable;

  Params()
      : timeLimit(900),
        treeSearch(1),
        onlyRelaxation(false),
        verbose(0),
        heuristicRootNode(3),
        heuristicOtherNodes(2),
        heuristic2stepVariant(3),
        heuristicSemigreedyAlpha(0.2),
        heuristicSemigreedyIter(100),
        feasibilityRootNode(2),
        feasibilityRootNodeTimeLimit(300),
        feasibilityOtherNodes(2),
        feasibilityOtherNodesTimeLimit(60),
        inheritColumns(0),
        initializationBigWeight(1000.0),
        preprocessing(true),
        pricingMethod(4),
        pricingHeur1Alpha(0.1),
        pricingHeur1MaxNCols(1),
        pricingMaxColsPerIter(10),
        pricingExactTimeLimit(300),
        branchingVariable(1) {};

  [[nodiscard]] bool is_verbose(int level = 1) const {
    return verbose >= level;
  }

  [[nodiscard]] bool use_dfs_tree_search() const { return treeSearch == 2; }

  [[nodiscard]] bool use_fms_branching() const {
    return branchingVariable == 1;
  }

  std::string get_pricing_method_name(int method) {
    switch (method) {
      case 0:
        return " (ILP)";
      case 1:
        return " (greedy -> ILP)";
      case 2:
        return " (greedy -> P,Q-MWSSP -> ILP)";
      case 3:
        return " (greedy -> P-MWSSP -> ILP)";
      case 4:
        return " (greedy -> P,Q-MWSSP -> P-MWSSP -> ILP)";
      case 5:
        return " (greedy -> P-MWSSP -> P,Q-MWSSP -> ILP)";
      case 6:
        return " (auto by density)";
      default:
        return " (unknown)";
    }
  };

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
      case 4:
        return " (auto)";
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

  std::string get_tree_search_name(int tree) {
    switch (tree) {
      case 1:
        return " (best-bound)";
      case 2:
        return " (dfs)";
      default:
        return " (unknown)";
    }
  };

  std::string get_branching_name(int branching) {
    switch (branching) {
      case 1:
        return " (fms)";
      case 2:
        return " (lntt)";
      default:
        return " (unknown)";
    }
  };

  void print_params(std::ostream& out) {
    out << "*** Parameters ***:" << std::endl;
    out << "Time limit: " << timeLimit << " seconds" << std::endl;
    out << "Tree search: " << treeSearch << get_tree_search_name(treeSearch)
        << std::endl;
    out << "Only relaxation: " << (onlyRelaxation ? "enabled" : "disabled")
        << std::endl;
    out << "Verbose level: " << verbose << std::endl;
    out << "Heuristic root node: " << heuristicRootNode
        << get_heur_name(heuristicRootNode);
    if (heuristicRootNode == 2 || heuristicRootNode == 3 ||
        heuristicRootNode == 4)
      out << ", variant: " << heuristic2stepVariant
          << get_heur_variant(heuristic2stepVariant);
    if (heuristicRootNode == 3 || heuristicRootNode == 4)
      out << ", alpha: " << heuristicSemigreedyAlpha
          << ", iterations: " << heuristicSemigreedyIter;
    out << std::endl;
    out << "Heuristic other nodes: " << heuristicOtherNodes
        << get_heur_name(heuristicOtherNodes);
    if (heuristicRootNode == 2 || heuristicRootNode == 3 ||
        heuristicRootNode == 4)
      out << ", variant: " << heuristic2stepVariant
          << get_heur_variant(heuristic2stepVariant);
    if (heuristicOtherNodes == 3 || heuristicRootNode == 4)
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
    out << "Preprocessing: " << (preprocessing ? "enabled" : "disabled")
        << std::endl;
    out << "Pricing method: " << pricingMethod
        << get_pricing_method_name(pricingMethod) << std::endl;
    out << "Branching variable: " << branchingVariable
        << get_branching_name(branchingVariable) << std::endl;
    out << "Pricing greedy max columns: " << pricingHeur1MaxNCols << std::endl;
    out << "Pricing max columns per iter: " << pricingMaxColsPerIter << std::endl;
    out << "Pricing heuristic 1 alpha: " << pricingHeur1Alpha << std::endl;
    out << "Pricing exact time limit: " << pricingExactTimeLimit << " seconds"
        << std::endl;
  }
};

#endif  // _PARAMS_HPP_