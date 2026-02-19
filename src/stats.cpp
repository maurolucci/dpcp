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
  size_t onodes = nodes == 1 ? 1 : nodes - 1;
  file << instance << "," << solver << "," << run << "," << nvertices << ","
       << nedges << "," << nA << "," << nB << "," << nvars << "," << ncons
       << "," << get_state_as_str() << "," << time << "," << nodes << ","
       << nodesLeft << "," << lb << "," << ub << "," << gap << ","
       << ninfeasPrepro + ninfeasCheck + ninfeasAux << "," << ninfeasPrepro
       << "," << ninfeasCheck << "," << ninfeasAux << "," << nint << "," << ngcp
       << "," << gcpTime / (ngcp > 0 ? ngcp : 1) << ","
       << nsolHeur + nsolLR + ntrivial << "," << nsolHeur << "," << nsolLR
       << "," << ntrivial << "," << bestTime << "," << bestIter << std::endl;

  file << rootlb << "," << rootub << "," << rootHeurTime << "," << rootFeasTime
       << ","
       << rootNCallsPool + rootNCallsHeur + rootNCallsMwis1 + rootNCallsMwis2 +
              rootNCallsExact
       << "," << rootNCallsPool << "," << rootNCallsHeur << ","
       << rootNCallsMwis1 << "," << rootNCallsMwis2 << "," << rootNCallsExact
       << ","
       << rootNColsPool + rootNColsHeur + rootNColsMwis1 + rootNColsMwis2 +
              rootNColsExact
       << "," << rootNColsPool << "," << rootNColsHeur << "," << rootNColsMwis1
       << "," << rootNColsMwis2 << "," << rootNColsExact << ","
       << rootTimePool + rootTimeHeur + rootTimeMwis1 + rootTimeMwis2 +
              rootTimeExact
       << "," << rootTimePool << "," << rootTimeHeur << "," << rootTimeMwis1
       << "," << rootTimeMwis2 << "," << rootTimeExact << std::endl;

  file << otherNodesHeurTime / onodes << "," << otherNodesFeasNCalls << ","
       << otherNodesFeasTime /
              (otherNodesFeasNCalls > 0 ? otherNodesFeasNCalls : 1)
       << ","
       << (otherNodesNCallsPool + otherNodesNCallsHeur + otherNodesNCallsMwis1 +
           otherNodesNCallsMwis2 + otherNodesNCallsExact) /
              onodes
       << "," << otherNodesNCallsPool / onodes << ","
       << otherNodesNCallsHeur / onodes << "," << otherNodesNCallsMwis1 / onodes
       << "," << otherNodesNCallsMwis2 / onodes << ","
       << otherNodesNCallsExact / onodes << ","
       << (otherNodesNColsPool + otherNodesNColsHeur + otherNodesNColsMwis1 +
           otherNodesNColsMwis2 + otherNodesNColsExact) /
              onodes
       << "," << otherNodesNColsPool / onodes << ","
       << otherNodesNColsHeur / onodes << "," << otherNodesNColsMwis1 / onodes
       << "," << otherNodesNColsMwis2 / onodes << ","
       << otherNodesNColsExact / onodes << ","
       << (otherNodesTimePool + otherNodesTimeHeur + otherNodesTimeMwis1 +
           otherNodesTimeMwis2 + otherNodesTimeExact) /
              onodes
       << "," << otherNodesTimePool / onodes << ","
       << otherNodesTimeHeur / onodes << "," << otherNodesTimeMwis1 / onodes
       << "," << otherNodesTimeMwis2 / onodes << ","
       << otherNodesTimeExact / onodes << std::endl;
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
  file << "Variables: " << nvars << std::endl;
  file << "Constraints: " << ncons << std::endl;
  file << "State: " << get_state_as_str() << std::endl;
  file << "Time total: " << time << std::endl;
  file << "Nodes total: " << nodes << std::endl;
  file << "Nodes left: " << nodesLeft << std::endl;
  file << "Lower bound: " << lb << std::endl;
  file << "Upper bound: " << ub << std::endl;
  file << "Gap: " << gap << std::endl;
  file << "Infeasible nodes: " << ninfeasPrepro + ninfeasCheck + ninfeasAux
       << " (total), " << ninfeasPrepro << " (preprocessing), " << ninfeasCheck
       << " (heuristic), " << ninfeasAux << " (auxiliary)" << std::endl;
  file << "Integer nodes: " << nint << std::endl;
  file << "GCP nodes: " << ngcp << std::endl;
  file << "GCP avg time: " << gcpTime / (ngcp > 0 ? ngcp : 1) << std::endl;
  file << "Solutions found: " << nsolHeur + nsolLR + ntrivial << " (total), "
       << nsolHeur << " (heuristic), " << nsolLR << " (linear relaxation), "
       << ntrivial << " (trivial)" << std::endl;
  file << "Best time: " << bestTime << std::endl;
  file << "Best iteration: " << bestIter << std::endl;
  file << "Root node stats:" << std::endl;
  file << "\tLower bound: " << rootlb << std::endl;
  file << "\tUpper bound: " << rootub << std::endl;
  file << "\tDPCP heuristic time: " << rootHeurTime << std::endl;
  file << "\tDPCP feasibility check time: " << rootFeasTime << std::endl;
  file << "\tPricing:" << std::endl;
  file << "\t\tCalls: "
       << rootNCallsPool + rootNCallsHeur + rootNCallsMwis1 + rootNCallsMwis2 +
              rootNCallsExact
       << " (total), " << rootNCallsPool << " (pool), " << rootNCallsHeur
       << " (heuristic), " << rootNCallsMwis1 << " (MWIS1), " << rootNCallsMwis2
       << " (MWIS2), " << rootNCallsExact << " (exact)" << std::endl;
  file << "\t\tColumns: "
       << rootNColsPool + rootNColsHeur + rootNColsMwis1 + rootNColsMwis2 +
              rootNColsExact
       << " (total), " << rootNColsPool << " (pool), " << rootNColsHeur
       << " (heuristic), " << rootNColsMwis1 << " (MWIS1), " << rootNColsMwis2
       << " (MWIS2), " << rootNColsExact << " (exact)" << std::endl;
  file << "\t\tTime: "
       << rootTimePool + rootTimeHeur + rootTimeMwis1 + rootTimeMwis2 +
              rootTimeExact
       << " (total), " << rootTimePool << " (pool), " << rootTimeHeur
       << " (heuristic), " << rootTimeMwis1 << " (MWIS1), " << rootTimeMwis2
       << " (MWIS2), " << rootTimeExact << " (exact)" << std::endl;
  size_t onodes = nodes == 1 ? 1 : nodes - 1;
  file << "Other nodes stats:" << std::endl;
  file << "\tDPCP heuristic avg time: " << otherNodesHeurTime / onodes
       << std::endl;
  file << "\tDPCP feasibility check calls: " << otherNodesFeasNCalls
       << std::endl;
  file << "\tDPCP feasibility check avg time: "
       << otherNodesFeasTime /
              (otherNodesFeasNCalls > 0 ? otherNodesFeasNCalls : 1)
       << std::endl;
  file << "\tPricing:" << std::endl;
  file << "\t\tCalls avg: "
       << (otherNodesNCallsPool + otherNodesNCallsHeur + otherNodesNCallsMwis1 +
           otherNodesNCallsMwis2 + otherNodesNCallsExact) /
              onodes
       << " (total), " << otherNodesNCallsPool / onodes << " (pool), "
       << otherNodesNCallsHeur / onodes << " (heuristic), "
       << otherNodesNCallsMwis1 / onodes << " (MWIS1), "
       << otherNodesNCallsMwis2 / onodes << " (MWIS2), "
       << otherNodesNCallsExact / onodes << " (exact)" << std::endl;
  file << "\t\tColumns avg: "
       << (otherNodesNColsPool + otherNodesNColsHeur + otherNodesNColsMwis1 +
           otherNodesNColsMwis2 + otherNodesNColsExact) /
              onodes
       << " (total), " << otherNodesNColsPool / onodes << " (pool), "
       << otherNodesNColsHeur / onodes << " (heuristic), "
       << otherNodesNColsMwis1 / onodes << " (MWIS1), "
       << otherNodesNColsMwis2 / onodes << " (MWIS2), "
       << otherNodesNColsExact / onodes << " (exact)" << std::endl;
  file << "\t\tTime avg: "
       << (otherNodesTimePool + otherNodesTimeHeur + otherNodesTimeMwis1 +
           otherNodesTimeMwis2 + otherNodesTimeExact) /
              onodes
       << " (total), " << otherNodesTimePool / onodes << " (pool), "
       << otherNodesTimeHeur / onodes << " (heuristic), "
       << otherNodesTimeMwis1 / onodes << " (MWIS1), "
       << otherNodesTimeMwis2 / onodes << " (MWIS2), "
       << otherNodesTimeExact / onodes << " (exact)" << std::endl;
}