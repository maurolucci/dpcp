#include "col.hpp"
#include "graph.hpp"
#include "stats.hpp"

using VertexSelector =
    std::function<Vertex(const GraphEnv &, const std::vector<bool> &)>;
using StableSetConstructor = std::function<VertexVector(
    const GraphEnv &, const std::vector<bool> &, const Vertex)>;

void dpcp_heur_1_step(const GraphEnv &genv, Col &col, VertexSelector getVertex,
                      StableSetConstructor buildStab);

Stats heur_solve(const GraphEnv &genv, const std::vector<TypeA> &as, Col &col,
                 size_t repetitions, Pool &pool);