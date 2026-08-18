#include "col.hpp"
#include "graph.hpp"
#include "params.hpp"
#include "stats.hpp"

Stats solve_ilp_cplex(DPCPInst& dpcp, const Params& params, std::ostream& log,
                      std::ostream& debugLog, Col& col);

Stats solve_ilp_gurobi(DPCPInst& dpcp, const Params& params, std::ostream& log,
                       std::ostream& debugLog, Col& col);