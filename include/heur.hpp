#include "col.hpp"
#include "graph.hpp"
#include "params.hpp"
#include "stats.hpp"

#include <iostream>

void dpcp_dsatur_heur(const GraphEnv &genv, VertexVector &selected,
                      std::map<TypeB, std::set<TypeB>> &adj, Col &col);

HeurStats dpcp_2_step_greedy_heur(const GraphEnv &genv, Col &col,
                                  const Params &params);

HeurStats dpcp_2_step_semigreedy_heur(const GraphEnv &genv, Col &col,
                                      const Params &params,
                                      std::ostream &iterFile);

HeurStats dpcp_1_step_greedy_heur(const GraphEnv &genv, Col &col);