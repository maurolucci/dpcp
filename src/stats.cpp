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
  case TIME_EXCEEDED_LP:
    return "TIME_EXCEEDED_LP";
  case TIME_EXCEEDED_PR:
    return "TIME_EXCEEDED_PR";
  case MEM_EXCEEDED:
    return "MEM_EXCEEDED";
  case MEM_EXCEEDED_LP:
    return "MEM_EXCEEDED_LP";
  case MEM_EXCEEDED_PR:
    return "MEM_EXCEEDED_PR";
  case INIT_FAIL:
    return "INIT_FAIL";
  default:
    return "UNKNOWN";
  }
}

void Stats::write_stats(std::ostream &file) {
  file << instance << "," << solver << "," << run << "," << nvertices << ","
       << nedges << "," << nA << "," << nB << "," << nvars << "," << ncons
       << "," << get_state_as_str() << "," << time << "," << nodes << ","
       << nodesLeft << "," << initSol << "," << initSolTime << "," << rootval
       << "," << lb << "," << ub << "," << gap << "," << ninfeas << "," << ngcp
       << "," << poolSize << "," << nCallsHeur << "," << nCallsMWis1 << ","
       << nCallsMWis2 << "," << nCallsExact << "," << nColsPool << ","
       << nColsHeur << "," << nColsMwis1 << "," << nColsMwis2 << ","
       << nColsExact << "," << nTimePool << "," << nTimeHeur << ","
       << nTimeMwis1 << "," << nTimeMwis2 << "," << nTimeExact << std::endl;
}

void Stats::print_stats(std::ostream &file) {
  file << std::endl << "*** Stats ***" << std::endl;
  file << "Instance: " << instance << std::endl;
  file << "Solver: " << solver << std::endl;
  file << "Run: " << run << std::endl;
  file << "Vertices: " << nvertices << std::endl;
  file << "Edges: " << nedges << std::endl;
  file << "|A|: " << nA << std::endl;
  file << "|B|: " << nB << std::endl;
  if (nvars != -1)
    file << "Variables: " << nvars << std::endl;
  if (ncons != -1)
    file << "Constraints: " << ncons << std::endl;
  file << "State: " << get_state_as_str() << std::endl;
  file << "Time: " << time << std::endl;
  file << "Nodes: " << nodes << std::endl;
  file << "Nodes left: " << nodesLeft << std::endl;
  file << "Initial solution: " << initSol << std::endl;
  file << "Initial solution time: " << initSolTime << std::endl;
  file << "Root relaxation: " << rootval << std::endl;
  file << "Lower bound: " << lb << std::endl;
  file << "Upper bound: " << ub << std::endl;
  file << "Gap: " << gap << std::endl;
  file << "Infeasible instances detected: " << ninfeas << std::endl;
  file << "GCP instances reached: " << ngcp << std::endl;
  file << "Size of pool: " << poolSize << std::endl;
  file << "Pricing (calls, cols, time):" << std::endl;
  file << "\tPool: " << nCallsPool << ", " << nColsPool << ", " << nTimePool
       << std::endl;
  file << "\tGreedy Heuristic: " << nCallsHeur << ", " << nColsHeur << ", "
       << nTimeHeur << std::endl;
  file << "\tMWSSP I: " << nCallsMWis1 << ", " << nColsMwis1 << ", "
       << nTimeMwis1 << std::endl;
  file << "\tMWSSP II: " << nCallsMWis2 << ", " << nColsMwis2 << ", "
       << nTimeMwis2 << std::endl;
  file << "\tExact: " << nCallsExact << ", " << nColsExact << ", " << nTimeExact
       << std::endl;
  file << "\tTotal time: "
       << nTimePool + nTimeHeur + nTimeMwis1 + nTimeMwis2 + nTimeExact
       << std::endl;
}