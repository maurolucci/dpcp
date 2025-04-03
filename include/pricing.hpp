#ifndef __PRICING_HPP__
#define __PRICING_HPP__

#include "graph.hpp"
#include "stats.hpp"

#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>

std::tuple<StableEnv, PRICING_STATE> exact_solve(GraphEnv &in,
                                                 IloNumArray &duals);

#endif // __PRICING_HPP__