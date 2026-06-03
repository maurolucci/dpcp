#include <iostream>
#include <limits>

#include "col.hpp"
#include "graph.hpp"
#include "params.hpp"
#include "stats.hpp"

void dpcp_dsatur_heur(const DPCPInst& dpcp, VertexVector& selected,
                      std::map<size_t, std::set<size_t>>& adj, Col& col);

HeurStats dpcp_2_step_greedy_heur(const DPCPInst& dpcp, Col& col,
                                  const Params& params);

HeurStats dpcp_2_step_semigreedy_heur(const DPCPInst& dpcp, Col& col,
                                      const Params& params,
                                      std::ostream& iterFile,
                                      double timelimit = 60.0);

HeurStats dpcp_2_step_semigreedy_heur(const DPCPInst& dpcp, Col& col,
                                      const Params& params,
                                      double timelimit = 60.0);

HeurStats dpcp_1_step_greedy_heur(const DPCPInst& dpcp, Col& col);