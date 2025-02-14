#include "lp.hpp"

#include <boost/graph/copy.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <cfloat>

LP::LP(const Graph &graph) : LP(Graph{graph}){};

LP::LP(const Graph &&graph) : graph(graph), mapA(), mapB(){};

// LP::LP(const HyperGraph &hg, const ConflictGraph &cg, double start_t)
//     : hg(hg), cg(cg), ecount(static_cast<int>(num_edges(cg))), elist(NULL),
//       nrows(hg.nbHyperedges() + hg.nbVertices()), start_t(start_t),
//       Xenv(IloEnv()), Xmodel(Xenv), Xvars(Xenv), Xrestr(Xenv),
//       Xobj(Xenv) /*, vars(Xenv)*/ {
//   vcount = num_vertices(cg);
//   vlist = new CGVertex[vcount];
//   size_t index = 0;
//   for (auto v : boost::make_iterator_range(vertices(cg)))
//     vlist[index++] = v;
//   index = 0;
//   elist = new int[2 * ecount];
//   auto [it, end] = edges(cg);
//   for (; it != end; ++it) {
//     elist[index++] = source(*it, cg);
//     elist[index++] = target(*it, cg);
//   }
// };

// LP::LP(const HyperGraph &hg, const ConflictGraph &cg, int vcount,
//        CGVertex *vlist, int ecount, int *elist, double start_t)
//     : hg(hg), cg(cg), ecount(ecount), elist(elist),
//       nrows(hg.nbHyperedges() + hg.nbVertices()), vcount(vcount),
//       vlist(vlist), start_t(start_t), Xenv(IloEnv()), Xmodel(Xenv),
//       Xvars(Xenv), Xrestr(Xenv), Xobj(Xenv) /*, vars(Xenv)*/ {};

// LP::~LP() {
//   // vars.end();
//   Xrestr.end();
//   Xvars.end();
//   Xobj.end();
//   Xmodel.end();
//   Xenv.end();
//   if (vlist != NULL)
//     delete[] vlist;
//   delete[] elist;
// }

void LP::initialize(CplexEnv &cenv) {

  // Fill mapA and mapB
  size_t iA = 0, iB = 0; // indices
  for (auto v : boost::make_iterator_range(vertices(graph))) {
    auto [a, b] = graph[v];
    if (mapA.insert({a, iA}).second)
      iA++;
    if (mapB.insert({b, iB}).second)
      iB++;
  }

  // Add constraints
  // Add ">= 1" constraints, one for each a \in A
  while (iA--)
    cenv.Xrestr.add(IloRange(cenv.Xenv, 1.0, IloInfinity));
  // Add "<= 1" constraints, one for each b \in B
  while (iB--)
    cenv.Xrestr.add(IloRange(cenv.Xenv, -1.0, IloInfinity));
  cenv.Xmodel.add(cenv.Xrestr);

  // Add objective function
  cenv.Xobj = IloMinimize(cenv.Xenv, 0.0);
  cenv.Xmodel.add(cenv.Xobj);

  // Add initial columns
  // Color each b \in B with an unique color
  for (auto &p : mapB) {
    TypeB b = p.first;
    std::set<TypeB> cB{b};
    // Find (*,b) vertices in graph
    std::set<TypeA> cA;
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      auto [a1, b1] = graph[v];
      if (b1 == b)
        cA.insert(a1);
    }
    add_column(cenv, cA, cB);
  }

  return;
}

void LP::add_column(CplexEnv &cenv, const std::set<TypeA> &cA,
                    const std::set<TypeB> &cB) {
  IloNumColumn column = cenv.Xobj(1.0);
  std::cout << "Adding column: [";
  for (auto a : cA) {
    column += cenv.Xrestr[mapA[a]](1.0);
    std::cout << " " << a;
  }
  std::cout << " ] [";
  for (auto b : cB) {
    column += cenv.Xrestr[mapA.size() + mapB[b]](-1.0);
    std::cout << " " << b;
  }
  std::cout << " ]" << std::endl;
  cenv.Xvars.add(IloNumVar(column));
}

void LP::add_column(CplexEnv &cenv, int count, const int *members) {
  std::set<TypeA> cA;
  std::set<TypeB> cB;
  for (int i = 0; i < count; ++i) {
    auto [a, b] = graph[members[i]];
    cA.insert(a);
    cB.insert(b);
  }
  add_column(cenv, cA, cB);
}

void LP::set_parameters(CplexEnv &cenv, IloCplex &cplex) {
  cplex.setDefaults();
  cplex.setParam(IloCplex::Param::Threads, 1);
  cplex.setParam(IloCplex::Param::Parallel, 1);
  cplex.setOut(cenv.Xenv.getNullStream());
  cplex.setWarning(cenv.Xenv.getNullStream());
}

LP_STATE LP::optimize() {

  auto startTime = std::chrono::high_resolution_clock::now();
  int rval = 0;
  LP_STATE state;

  CplexEnv cenv;

  MWISenv *mwis_env = NULL;
  COLORNWT *mwis_pi = NULL;
  COLORNWT mwis_pi_scalef = 1;

  COLORset *newsets = NULL;
  int nnewsets = 0;

  // TODO: Check
  // If some hyperedge has no representator -> INFLEASIBLE
  // If every hyperedge has exactly one representator -> SOLVE COLORING
  // Otherwise, do the code below

  // Initialize cplex environment
  IloCplex cplex(cenv.Xmodel);
  set_parameters(cenv, cplex);

  // Initialize linear relaxation
  initialize(cenv);

  // Initalize stable environment
  rval = COLORstable_initenv(&mwis_env, NULL, 0);
  mwis_pi = (COLORNWT *)COLOR_SAFE_MALLOC(num_vertices(graph), COLORNWT);

  do {

    // Reset solution counter
    nnewsets = 0;

    // Set time limit
    auto currentTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration(currentTime - startTime);
    double timeLimit = TIMELIMIT - duration.count();
    if (timeLimit < 0) {
      state = TIME_OR_MEM_LIMIT;
      break;
    }
    cplex.setParam(IloCplex::Param::TimeLimit, timeLimit);

    // Optimize
    cplex.solve();

    // Handle errores
    IloCplex::CplexStatus status = cplex.getCplexStatus();
    if (status == IloCplex::AbortTimeLim || status == IloCplex::MemLimFeas ||
        status == IloCplex::MemLimInfeas) {
      state = TIME_OR_MEM_LIMIT;
      break;
    }

    // Recover dual values
    IloNumArray duals(cenv.Xenv, mapA.size() + mapB.size());
    cplex.getDuals(duals, cenv.Xrestr);
    std::cout << "duals: ";
    for (int i = 0; i < mapA.size() + mapB.size(); ++i)
      std::cout << duals[i] << " ";
    std::cout << std::endl;
    IloNumArray weights(cenv.Xenv, num_vertices(graph));
    size_t i = 0;
    for (auto v : boost::make_iterator_range(vertices(graph))) {
      auto [a, b] = graph[v];
      weights[i++] = duals[mapA[a]] - duals[mapA.size() + mapB[b]];
    }

    // Scale weights
    std::cout << "weights: ";
    for (int i = 0; i < num_vertices(graph); ++i)
      std::cout << weights[i] << " ";
    std::cout << std::endl;
    double2COLORNWT(mwis_pi, &mwis_pi_scalef, weights);

    // Get edge list
    int ecount = 0;
    int elist[2 * num_edges(graph)];
    for (auto e : boost::make_iterator_range(edges(graph))) {
      elist[2 * ecount] = source(e, graph);
      elist[2 * ecount++ + 1] = target(e, graph);
    }

    // Solve the MWIS problem
    rval =
        COLORstable_wrapper(&mwis_env, &newsets, &nnewsets, num_vertices(graph),
                            ecount, elist, mwis_pi, mwis_pi_scalef, 0, 0, 2);

    // Add column/s
    for (int set_i = 0; set_i < nnewsets; ++set_i) {
      add_column(cenv, newsets[set_i].count, newsets[set_i].members);
      free(newsets[set_i].members);
    }
    if (newsets != NULL)
      free(newsets);

  } while (nnewsets > 0);

  if (state != TIME_OR_MEM_LIMIT) {

    // Recover primal values and objective value
    IloNumArray values = IloNumArray(cenv.Xenv, cenv.Xvars.getSize());
    cplex.getValues(values, cenv.Xvars);
    std::cout << "Primal values: ";
    for (int i = 0; i < cenv.Xvars.getSize(); ++i)
      std::cout << values[i] << " ";
    std::cout << std::endl;
    IloNum obj_value = cplex.getObjValue();
    std::cout << "LR value: " << obj_value << std::endl;

    // Check for integrality
    bool is_integer = true;
    for (int i = 0; i < cenv.Xvars.getSize(); ++i)
      if (values[i] < 1 - EPSILON)
        is_integer = false;

    if (is_integer) {
      obj_value = round(obj_value);
      state = INTEGER;
    }
    state = FRACTIONAL;

    values.end();
  }

  free(mwis_pi);
  COLORstable_freeenv(&mwis_env);
  cplex.end();
  return state;
}

/* Adaptation of  COLOR_double2COLORNWT from exactcolors.
The type of dbl_nweights change from const double[] to const IloNumArray.
*/
int LP::double2COLORNWT(COLORNWT nweights[], COLORNWT *scalef,
                        const IloNumArray dbl_nweights) {
  size_t i;
  double max_dbl_nweight = -DBL_MAX;
  double max_prec_dbl = exp2(DBL_MANT_DIG - 1);
  static const double max_mwiswt = (double)COLORNWT_MAX;
  double dbl_scalef = COLORDBLmin(max_prec_dbl, max_mwiswt);

  dbl_scalef /= (double)dbl_nweights.getSize();

  for (i = 0; i < dbl_nweights.getSize(); ++i) {
    max_dbl_nweight = COLORDBLmax(max_dbl_nweight, dbl_nweights[i]);
  }
  dbl_scalef /= COLORDBLmax(1.0, max_dbl_nweight);
  dbl_scalef = floor(dbl_scalef);
  *scalef = (COLORNWT)dbl_scalef;

  for (i = 0; i < dbl_nweights.getSize(); ++i) {
    double weight = dbl_nweights[i] * dbl_scalef;
    assert(weight < (double)COLORNWT_MAX);
    nweights[i] = (COLORNWT)weight;
  }
  return 0;
}

void LP::save_solution(Coloring &col) {
  /*
    value = G->get_precoloring_value();
    std::vector<int> stables_per_color(G->get_n_colors(), 0);
    std::vector<int> temp_coloring(G->get_n_vertices());

    // Build the coloring of the current graph
    for (int i : pos_vars) {
      int color = vars[i].color;
      int true_color = G->get_C(color, stables_per_color[color]);
      stables_per_color[color]++;
      value += G->get_color_cost(color);
      for (int j = 0; j < G->get_n_vertices(); ++j)
        if (vars[i].stable[j])
          temp_coloring[j] = true_color;
    }

    // Build the coloring of the original graph
    coloring.resize(G->get_n_total_vertices());
    for (int i = 0; i < G->get_n_total_vertices(); ++i) {
      int cv = G->get_current_vertex(i);
      if (cv == -1)
        coloring[i] = G->get_precoloring(i);
      else
        coloring[i] = temp_coloring[cv];
      active_colors.insert(coloring[i]);
    }
  */
  return;
}

void LP::branch(std::vector<LP *> &branches, Vertex v) {

  if (N_BRANCHES != 2) {
    std::cout << "N_BRANCHES != 2: Unimplemented" << std::endl;
    abort();
  }

  // Recover a and b
  auto [a, b] = graph[v];

  // *******
  // ** Left branch: (a,b) is colored
  // *******

  // Get the set of vertices to be removed
  // I.e. every vertex (a,b') with b' != b
  std::set<Vertex> removed;
  for (auto v : boost::make_iterator_range(vertices(graph)))
    if (graph[v].first == a && graph[v].second != b)
      removed.insert(v);

  // Create a view of the graph without the vertices to be removed
  using VFilter = boost::is_not_in_subset<std::set<Vertex>>;
  VFilter vFilter(removed);
  using Filter = boost::filtered_graph<Graph, boost::keep_all, VFilter>;
  Filter filtered(graph, boost::keep_all(), vFilter);

  // Create a copy of the view (force reindex)
  Graph graph1;
  boost::copy_graph(filtered, graph1);

  // *******
  // ** Right branch: (a,b) is uncolored
  // ** ¡Reuse graph!
  // *******

  Graph graph2(std::move(graph));
  remove_vertex(v, graph2);

  // *******
  // ** Create branches
  // *******

  branches.resize(2);
  branches[0] = new LP(graph1);
  branches[1] = new LP(graph2);

  for (auto v : boost::make_iterator_range(vertices(graph1)))
    std::cout << v << ": (" << graph1[v].first << "," << graph1[v].second
              << ") ";
  std::cout << std::endl;

  for (auto v : boost::make_iterator_range(vertices(graph2)))
    std::cout << v << ": (" << graph2[v].first << "," << graph2[v].second
              << ") ";
  std::cout << std::endl;

  return;
}

/*
void LP::branch(std::vector<LP *> &branches, Vertex v) {

  if (N_BRANCHES != 2) {
    std::cout << "N_BRANCHES != 2: Unimplemented" << std::endl;
    abort();
  }

  // Recover a and b
  auto [a, b] = graph[v];

  // Find some index of interest
  // .. (e-1,u), (e,v1), ..., (e,v), ..., (e,vn), (e+1,w)
  //               i1           i2                  i3
  size_t i1 = 0, i2 = 0, i3 = 0;
  for (; i1 < vcount; ++i1)
    if (cg[vlist[i1]].first == e)
      break;
  for (i2 = i1; i2 < vcount; ++i2)
    if (cg[vlist[i2]].second == v)
      break;
  for (i3 = i2 + 1; i3 < vcount; ++i3)
    if (cg[vlist[i3]].first != e)
      break;

  // *******
  // ** Left branch: e is represented with v
  // *******

  int vcount1 = vcount - i3 + i1 + 1;
  CGVertex *vlist1 = new CGVertex[2 * vcount1];
  memcpy(vlist1, vlist, sizeof(CGVertex) * i1);
  vlist1[i1] = cgv;
  memcpy(vlist1 + i1 + 1, vlist + i3, sizeof(CGVertex) * (vcount - i3));

  int ecount1 = 0;
  int *elist1 = new int[2 * ecount];
  for (int i = 0; i < 2 * ecount; ++i) {
    int u = elist[i];
    if (u >= i1 && u < i3 && u != i2) {
      if (i % 2 == 0)
        ++i;
      continue;
    }
    int u2 = u < i1 ? u : (u == i2 ? i1 : u - i3 + i1 + 1);
    if (i % 2 == 0)
      elist1[2 * ecount1] = u2;
    else {
      elist1[2 * ecount1 + 1] = u2;
      ecount1++;
    }
  }

  LP *lp1 = new LP(hg, cg, vcount1, vlist1, ecount1, elist1, start_t);

  // *******
  // ** Right branch: e is not represented with v
  // *******

  int vcount2 = vcount - 1;
  CGVertex *vlist2 = vlist;
  vlist = NULL;
  memcpy(vlist2 + i2, vlist2 + i2 + 1, sizeof(CGVertex) * (vcount - i2 - 1));

  int ecount2 = 0;
  int *elist2 = new int[2 * ecount];
  for (int i = 0; i < 2 * ecount; ++i) {
    int u = elist[i];
    if (u == i2) {
      if (i % 2 == 0)
        ++i;
      continue;
    }
    int u2 = u < i2 ? u : u - 1;
    if (i % 2 == 0)
      elist2[2 * ecount2] = u2;
    else {
      elist2[2 * ecount2 + 1] = u2;
      ecount2++;
    }
  }

  LP *lp2 = new LP(hg, cg, vcount2, vlist2, ecount2, elist2, start_t);

  branches.resize(2);
  branches[0] = lp1;
  branches[1] = lp2;

  // TODO: Remove
  std::cout << "i1, i2, i3: " << i1 << " " << i2 << " " << i3 << std::endl;

  std::cout << "vcount1: " << vcount1 << std::endl;
  std::cout << "vlist1: ";
  for (int i = 0; i < vcount1; ++i)
    std::cout << i << ":(" << cg[vlist1[i]].first << "," <<
  cg[vlist1[i]].second
              << ") ";
  std::cout << std::endl;

  std::cout << "ecount1: " << ecount1 << std::endl;
  std::cout << "elist1: ";
  for (int i = 0; i < ecount1; ++i)
    std::cout << "(" << elist1[2 * i] << "," << elist1[2 * i + 1] << ") ";
  std::cout << std::endl;

  std::cout << "vcount2: " << vcount2 << std::endl;
  std::cout << "vlist2: ";
  for (int i = 0; i < vcount2; ++i)
    std::cout << i << ":(" << cg[vlist2[i]].first << "," <<
  cg[vlist2[i]].second
              << ") ";
  std::cout << std::endl;

  std::cout << "ecount2: " << ecount2 << std::endl;
  std::cout << "elist2: ";
  for (int i = 0; i < ecount2; ++i)
    std::cout << "(" << elist2[2 * i] << "," << elist2[2 * i + 1] << ") ";
  std::cout << std::endl;

return;
}
*/