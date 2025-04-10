#include "col.hpp"
#include "graph.hpp"
#include "stats.hpp"

Stats heur_solve(const GraphEnv &genv, const std::vector<TypeA> &as, Col &col,
                 size_t repetitions);