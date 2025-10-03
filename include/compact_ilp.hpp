#include "col.hpp"
#include "graph.hpp"
#include "params.hpp"
#include "stats.hpp"

Stats solve_ilp(const Graph &graph, const Params &params, std::ostream &log,
                Col &col);