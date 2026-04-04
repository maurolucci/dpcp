#include "bp.hpp"

#include <cmath>
#include <iomanip>
#include <optional>

namespace {
constexpr double kGapDenominatorEpsilon = 1e-6;
constexpr double kLogIntervalSeconds = 10.0;

inline bool should_prune_by_bound(double lowerBound, double primalBound) {
  return std::ceil(lowerBound - EPSILON_BP) >= primalBound;
}
}  // namespace

Node::Node(std::unique_ptr<LP> lp, size_t depth, size_t id)
    : lp(std::move(lp)), depth(depth), id(id) {}

Node::Node(LP&& lp, size_t depth, size_t id)
    : Node(std::make_unique<LP>(std::move(lp)), depth, id) {}

double Node::get_obj_value() const { return lp->get_lower_bound(); }

bool Node::operator>(const Node& n) const {
  return (get_obj_value() > n.get_obj_value());
}

LP_STATE Node::solve(double timelimit, double ub) {
  return lp->solve(timelimit, ub);
}

bool Node::feas_sol() const { return lp->has_heur_solution(); }

size_t Node::feas_value() const { return lp->get_upper_bound(); }

LP_INTEGER_SOURCE Node::integer_source() const {
  return lp->get_integer_source();
}

void Node::save(Col& sol) { sol = lp->get_lp_solution(); }

void Node::save_heur(Col& sol) { sol = lp->get_heur_solution(); }

void Node::branch(std::vector<Node>& sons) {
  std::vector<std::unique_ptr<LP>> childLps;
  lp->branch(childLps);
  sons.clear();
  sons.reserve(childLps.size());
  for (auto& childLp : childLps) {
    sons.emplace_back(std::move(childLp), depth + 1);
  }
}

BP::BP(Params& params, std::ostream& log, std::ostream& debugLog, Col& sol,
       double ub)
    : params(params),
      best_integer_solution(sol),
      primal_bound(ub),
      nodes(0),
      nextNodeId(0),
      first_call(true),
      log(log),
      debugLog(debugLog),
      stats() {}

Stats BP::solve(DPCPInst& originalDpcp) {
  start_t = ClockType::now();
  last_t = start_t;
  first_call = true;
  nextNodeId = 0;

  // Make a copy of the original DPCP instance to preprocess and modify during
  // the algorithm, while keeping the original one unchanged for reference
  DPCPInst dpcp(originalDpcp);
  log << "Original DPCP instance: |V|=" << num_vertices(origDpcp.get_graph())
      << ", |E|=" << num_edges(origDpcp.get_graph())
      << ", |P|=" << origDpcp.get_nP() << ", |Q|=" << origDpcp.get_nQ()
      << std::endl;
  if (params.preprocessing) dpcp.preprocess(true);
  log << "After preprocessing: |V|=" << num_vertices(dpcp.get_graph())
      << ", |E|=" << num_edges(dpcp.get_graph()) << ", |P|=" << dpcp.get_nP()
      << ", |Q|=" << dpcp.get_nQ() << std::endl;
  Pool pool;
  LP rootLp(std::move(dpcp), std::move(pool), origDpcp, params, stats, log,
            debugLog, true);
  Node root(std::move(rootLp), 0, nextNodeId++);

  auto state_after_push = [](LP_STATE lpState) -> std::optional<STATE> {
    switch (lpState) {
      case LP_TIME_EXCEEDED:
        return TIME_EXCEEDED_LP;
      case LP_TIME_EXCEEDED_PR:
        return TIME_EXCEEDED_PR;
      case LP_MEM_EXCEEDED:
        return MEM_EXCEEDED_LP;
      case LP_MEM_EXCEEDED_PR:
        return MEM_EXCEEDED_PR;
      case LP_INIT_FAIL:
        return INIT_FAIL;
      default:
        return std::nullopt;
    }
  };

  // Push root node (and solve initial LR)
  if (const auto state = state_after_push(push(std::move(root)));
      state.has_value()) {
    return return_stats(*state);
  }

  if (!params.onlyRelaxation) {
    while (!L.empty()) {
      // Pop next node
      show_stats();  // First show_stats, then pop
      Node node = std::move(L.back());
      pop();

      // Re-try to prune by bound, since primal_bound could have been improved
      if (ceil(node.get_obj_value() - EPSILON_BP) >= primal_bound) continue;

      // Branch
      std::vector<Node> sons;
      node.branch(sons);

      if (params.is_verbose(2)) {
        debugLog << "Branch of node id=" << node.get_id()
                 << " at depth=" << node.get_depth() << " added " << sons.size()
                 << " sons" << std::endl;
      }

      // Push sons (and solve initial LR)
      for (auto& n : sons) {
        n.set_id(nextNodeId++);
        if (const auto state = state_after_push(push(std::move(n)));
            state.has_value()) {
          return return_stats(*state);
        }
      }
    }
  }

  if (primal_bound == DBL_MAX) return return_stats(INFEASIBLE);

  if (params.onlyRelaxation && !L.empty()) return return_stats(FEASIBLE);

  return return_stats(OPTIMAL);
}

size_t BP::get_nodes() { return nodes; }

double BP::get_gap() {
  const double dual_bound = calculate_dual_bound();
  double gap = DBL_MAX;
  if (dual_bound != -DBL_MAX)
    gap = std::fabs(dual_bound - primal_bound) /
          (kGapDenominatorEpsilon + std::fabs(primal_bound)) * 100;
  return gap;
}

double BP::get_primal_bound() { return primal_bound; }

Stats BP::return_stats(STATE state) {
  stats.state = state;
  stats.time =
      std::chrono::duration<double>(ClockType::now() - start_t).count();
  if (stats.time > params.timeLimit) {
    stats.time = params.timeLimit;
    stats.state = primal_bound == DBL_MAX ? TIME_EXCEEDED : FEASIBLE;
  }
  stats.nodes = static_cast<int>(nodes);
  stats.nodesLeft = static_cast<int>(L.size());
  stats.ub = static_cast<int>(primal_bound + 0.5);
  if (state == OPTIMAL) {
    stats.lb = primal_bound;
    stats.gap = 0.0;
  } else if (state != INFEASIBLE) {
    stats.lb = calculate_dual_bound();
    stats.gap = get_gap() / 100;
  }

  return stats;
}

LP_STATE BP::push(Node node) {
  nodes++;
  double obj_value;

  if (params.is_verbose(2)) {
    debugLog << "Solving node id=" << node.get_id()
             << " depth=" << node.get_depth() << std::endl;
  }

  // Solve the linear relaxation of the node and prune if possible
  const double elapsed =
      std::chrono::duration<double>(ClockType::now() - start_t).count();
  LP_STATE state = node.solve(params.timeLimit - elapsed, primal_bound);

  // If a better feasible solution was found, save it
  if (node.feas_sol()) {
    obj_value = node.feas_value();
    if (obj_value < primal_bound) {
      stats.nsolHeur++;
      node.save_heur(best_integer_solution);
      log << "New best integer solution found by heuristic with value: "
          << obj_value << std::endl;
      update_primal_bound(obj_value);
    }
  }

  switch (state) {
    case LP_INTEGER:
      obj_value = node.get_obj_value();
      if (obj_value < primal_bound) {
        node.save(best_integer_solution);
        switch (node.integer_source()) {
          case LP_INTEGER_SOURCE_GCP:
            stats.nsolGCP++;
            log << "New best integer solution found by GCP with value: "
                << obj_value << std::endl;
            break;
          case LP_INTEGER_SOURCE_TRIVIAL:
            stats.nsolTrivial++;
            log << "New best integer solution found by trivial reduction with "
                   "value: "
                << obj_value << std::endl;
            break;
          default:
            stats.nsolLR++;
            log << "New best integer solution found by LR with value: "
                << obj_value << std::endl;
            break;
        }
        update_primal_bound(obj_value);
      }
      return state;  // Prune by optimality

    case LP_FRACTIONAL:
      obj_value = node.get_obj_value();
      if (should_prune_by_bound(obj_value, primal_bound))
        return state;  // Prune by bound
      break;           // Do not prune

    default:
      return state;  // Prune by infeasibility or mem/time limit
  }

  if (params.use_dfs_tree_search()) {
    L.push_back(std::move(node));
    return LP_FRACTIONAL;
  }

  if (L.empty()) {
    L.push_back(std::move(node));
    return state;
  }

  const double node_obj_value = node.get_obj_value();
  for (auto it = L.begin(); it != L.end(); ++it)
    if (node_obj_value > it->get_obj_value()) {
      L.insert(it, std::move(node));
      return state;
    }

  L.push_back(std::move(node));
  return state;
}

Node& BP::top() { return L.back(); }

void BP::pop() { L.pop_back(); }

void BP::update_primal_bound(double obj_value) {
  primal_bound = obj_value;

  if (params.use_dfs_tree_search()) {
    for (auto it = L.begin(); it != L.end();) {
      const double node_obj_value = it->get_obj_value();
      if (node_obj_value >= primal_bound)
        it = L.erase(it);
      else
        ++it;
    }
    return;
  }

  for (auto it = L.begin(); it != L.end();) {
    const double node_obj_value = it->get_obj_value();
    if (node_obj_value >= primal_bound)
      it = L.erase(it);
    else
      break;  // Since the list is sorted, we can stop as soon as we find a node
              // that cannot be pruned
  }
}

double BP::calculate_dual_bound() {
  double dual_bound = DBL_MAX;
  if (!L.empty()) {
    if (params.use_dfs_tree_search()) {
      for (auto it = L.begin(); it != L.end(); ++it)
        if (it->get_obj_value() < dual_bound) dual_bound = it->get_obj_value();
    } else {
      dual_bound = L.front().get_obj_value();
    }
  }

  if (dual_bound == DBL_MAX) {
    return -DBL_MAX;
  }
  return dual_bound;
}

void BP::show_stats() {
  auto now_t = ClockType::now();

  if (first_call)
    first_call = false;
  else {
    if (std::chrono::duration<double>(now_t - last_t).count() <
        kLogIntervalSeconds)
      return;
    last_t = now_t;
  }

  const double dual_bound = calculate_dual_bound();

  log << std::fixed << std::setprecision(3);

  log << "  Best obj value  = " << dual_bound << "\t Best int = ";
  if (primal_bound == DBL_MAX)
    log << "inf\t Gap = ---";
  else
    log << (int)(EPSILON_BP + primal_bound) << "\t Gap = "
        << (primal_bound - dual_bound) / (EPSILON_BP + primal_bound) * 100
        << "%";
  log << "\t Nodes: processed = " << nodes << ", left = " << L.size()
      << "\t time = " << std::chrono::duration<double>(now_t - start_t).count()
      << std::endl;
}
