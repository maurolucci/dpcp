#include "col.hpp"
#include "graph.hpp"
#include "stats.hpp"

Stats solve_ilp(const Graph &graph, size_t ncolors, std::ostream &log,
                const Col &initialCol);