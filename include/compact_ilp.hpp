#include "col.hpp"
#include "graph.hpp"
#include "params.hpp"
#include "stats.hpp"

Stats solve_ilp(DPCPInst &dpcp, const Params &params, std::ostream &log,
                std::ostream &debugLog, Col &col);