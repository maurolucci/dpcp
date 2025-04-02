#include "stats.hpp"

std::string Stats::get_state_as_str() {
  switch (state) {
  case OPTIMAL:
    return "OPTIMAL";
  case FEASIBLE:
    return "FEASIBLE";
  case INFEASIBLE:
    return "INFEASIBLE";
  case TIME_EXCEEDED:
    return "TIME_EXCEEDED";
  case MEM_EXCEEDED:
    return "MEM_EXCEEDED";
  case NODE_TIME_EXCEEDED:
    return "NODE_TIME_EXCEEDED";
  case NODE_MEM_EXCEEDED:
    return "NODE_MEM_EXCEEDED";
  default:
    return "UNKNOWN";
  }
}

void Stats::write_stats(std::ostream &file) {
  file << nvars << "," << ncons << "," << get_state_as_str() << "," << time
       << "," << nodes << "," << lb << "," << ub << "," << gap << std::endl;
}

void Stats::print_stats(std::ostream &file) {
  file << std::endl << "*** Stats ***" << std::endl;
  file << "Variables: " << nvars << std::endl;
  file << "Constraints: " << ncons << std::endl;
  file << "State: " << get_state_as_str() << std::endl;
  file << "Time: " << time << std::endl;
  file << "Nodes: " << nodes << std::endl;
  file << "Lower bound: " << lb << std::endl;
  file << "Upper bound: " << ub << std::endl;
  file << "Gap: " << gap << std::endl;
}