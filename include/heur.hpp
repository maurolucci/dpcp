#include <iostream>
#include <limits>
#include <optional>
#include <unordered_map>
#include <unordered_set>

#include "col.hpp"
#include "graph.hpp"
#include "params.hpp"
#include "stats.hpp"

using Column = StableEnv;
using Pool = std::vector<Column>;

HeurStats dpcp_2_step_greedy_heur(const DPCPInst& dpcp, Col& col,
                                  const Params& params);

HeurStats dpcp_2_step_semigreedy_heur(const DPCPInst& dpcp, Col& col,
                                      const Params& params,
                                      std::ostream& iterFile,
                                      Pool* pool = nullptr);

HeurStats dpcp_2_step_semigreedy_heur(const DPCPInst& dpcp, Col& col,
                                      const Params& params,
                                      Pool* pool = nullptr);

HeurStats dpcp_1_step_greedy_heur(const DPCPInst& dpcp, Col& col);