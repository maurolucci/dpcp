#include "lp.hpp"

#include <cfloat>

// Column::Column(int n_best, const nodepnt *best_sol) {}

LP::LP(const HyperGraph &hg, const ConflictGraph &cg, double start_t)
    : hg(hg), cg(cg), ecount(static_cast<int>(num_edges(cg))), elist(NULL),
      nrows(hg.nbHyperedges() + hg.nbVertices()), start_t(start_t),
      Xenv(IloEnv()), Xmodel(Xenv), Xvars(Xenv), Xrestr(Xenv),
      Xobj(Xenv) /*, vars(Xenv)*/ {
  vcount = num_vertices(cg);
  vlist = new CGVertex[vcount];
  size_t index = 0;
  for (auto v : boost::make_iterator_range(vertices(cg)))
    vlist[index++] = v;
  index = 0;
  elist = new int[2 * ecount];
  auto [it, end] = edges(cg);
  for (; it != end; ++it) {
    elist[index++] = source(*it, cg);
    elist[index++] = target(*it, cg);
  }
};

LP::LP(const HyperGraph &hg, const ConflictGraph &cg, int vcount,
       CGVertex *vlist, int ecount, int *elist, double start_t)
    : hg(hg), cg(cg), ecount(ecount), elist(elist),
      nrows(hg.nbHyperedges() + hg.nbVertices()), vcount(vcount), vlist(vlist),
      start_t(start_t), Xenv(IloEnv()), Xmodel(Xenv), Xvars(Xenv), Xrestr(Xenv),
      Xobj(Xenv) /*, vars(Xenv)*/ {};

LP::~LP() {
  // vars.end();
  Xrestr.end();
  Xvars.end();
  Xobj.end();
  Xmodel.end();
  Xenv.end();
  if (vlist != NULL)
    delete[] vlist;
  delete[] elist;
}

void LP::initialize() {
  // Initialize hyperedge and vertex constraints
  // We will have "hyperedge" constraints with r.h.s >= 1 and "vertex"
  // constraints with r.h.s >= -1
  for (int i = 0; i < hg.nbHyperedges(); i++)
    Xrestr.add(IloRange(Xenv, 1.0, IloInfinity));
  for (int i = 0; i < hg.nbVertices(); i++)
    Xrestr.add(IloRange(Xenv, -1.0, IloInfinity));
  Xmodel.add(Xrestr);

  // Initialize objective function
  Xobj = IloMinimize(Xenv, 0.0);
  Xmodel.add(Xobj);

  fill_initial_columns();

  return;
}

void LP::fill_initial_columns() {
  // TODO: Mejorar este algoritmo
  for (const auto &v : hg.vertices()) {
    int count = 0;
    int members[vcount];
    for (int i = 0; i < vcount; ++i) {
      auto [e, u] = cg[vlist[i]];
      if (u == v->id())
        members[count++] = i;
    }
    if (count > 0)
      add_column(count, members);
  }
  return;
}

void LP::add_column(int count, const int *members) {
  IloNumColumn column = Xobj(1.0);
  std::vector<bool> hgCoef(hg.nbHyperedges(), false);
  std::vector<bool> vCoef(hg.nbVertices(), false);
  std::cout << "Adding column: ";
  for (int i = 0; i < count; ++i) {
    std::cout << members[i];
    auto [e, v] = cg[vlist[members[i]]];
    std::cout << "(" << e << "," << v << "), ";
    // Fill the column corresponding to ">= 1" constraints
    if (!hgCoef[e]) {
      column += Xrestr[e](1.0);
      hgCoef[e] = true;
    }
    // and the ">= -1 constraint
    if (!vCoef[v]) {
      column += Xrestr[hg.nbHyperedges() + v](-1.0);
      vCoef[v] = true;
    }
  }
  std::cout << std::endl;
  // add the column as a non-negative continuos variable
  Xvars.add(IloNumVar(column));
}

void LP::set_parameters(IloCplex &cplex) {
  cplex.setDefaults();
#ifndef SHOWCPLEX
  cplex.setOut(Xenv.getNullStream());
  cplex.setWarning(Xenv.getNullStream());
#endif
  cplex.setParam(IloCplex::Param::Threads, 1);
  cplex.setParam(IloCplex::Param::Parallel, 1);
}

LP_STATE LP::optimize() {

  auto start_time = std::chrono::high_resolution_clock::now();
  int rval = 0;

  MWISenv *mwis_env = NULL;
  COLORNWT *mwis_pi = NULL;
  COLORNWT mwis_pi_scalef = 1;

  COLORset *newsets = NULL;
  int nnewsets = 0;

  // TODO: Check
  // If some hyperedge has no representator -> INFLEASIBLE
  // If every hyperedge has exactly one representator -> SOLVE COLORING
  // Otherwise, do the code below

  IloCplex cplex(Xmodel);
  set_parameters(cplex);
  initialize();

  // Arguments:
  //  pname = NULL;
  //  write_mwis = 0;
  rval = COLORstable_initenv(&mwis_env, NULL, 0);

  mwis_pi = (COLORNWT *)COLOR_SAFE_MALLOC(vcount, COLORNWT);

  do {

    // Reset solution counter
    nnewsets = 0;

    // Set time limit
    auto current_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration(current_time - start_time);
    double time_limit = MAXTIME - duration.count();
    if (time_limit < 0) {
      free(mwis_pi);
      COLORstable_freeenv(&mwis_env);
      cplex.end();
      return TIME_OR_MEM_LIMIT;
    }
    cplex.setParam(IloCplex::Param::TimeLimit, time_limit);

    // Optimization
    cplex.solve();

    // Handle non-optimality
    IloCplex::CplexStatus status = cplex.getCplexStatus();
    if (status == IloCplex::AbortTimeLim || status == IloCplex::MemLimFeas ||
        status == IloCplex::MemLimInfeas) {
      free(mwis_pi);
      COLORstable_freeenv(&mwis_env);
      cplex.end();
      return TIME_OR_MEM_LIMIT;
    }

    // Recover dual values
    IloNumArray duals(Xenv, nrows);
    cplex.getDuals(duals, Xrestr);
    std::cout << "duals: ";
    for (int i = 0; i < nrows; ++i)
      std::cout << duals[i] << " ";
    std::cout << std::endl;
    IloNumArray weights(Xenv, vcount);
    for (size_t i = 0; i < vcount; ++i) {
      auto [e, v] = cg[vlist[i]];
      weights[i] = duals[e] - duals[hg.nbHyperedges() + v];
    }

    // Scale weights
    std::cout << "weights: ";
    for (int i = 0; i < vcount; ++i)
      std::cout << weights[i] << " ";
    std::cout << std::endl;
    double2COLORNWT(mwis_pi, &mwis_pi_scalef, weights);

    // Solve the MWIS problem. Arguments:
    //  greedy_only = 0
    //  force_rounding = 0
    //  rounding_strategy = 2 (no_rounding)
    rval = COLORstable_wrapper(&mwis_env, &newsets, &nnewsets, vcount, ecount,
                               elist, mwis_pi, mwis_pi_scalef, 0, 0, 2);

    // Add column/s
    for (int set_i = 0; set_i < nnewsets; ++set_i) {
      add_column(newsets[set_i].count, newsets[set_i].members);
      free(newsets[set_i].members);
    }
    if (newsets != NULL)
      free(newsets);

  } while (nnewsets > 0);

  // Recover primal values and objective value
  values = IloNumArray(Xenv, Xvars.getSize());
  cplex.getValues(values, Xvars);
  std::cout << "Primal values: ";
  for (int i = 0; i < Xvars.getSize(); ++i)
    std::cout << values[i] << " ";
  std::cout << std::endl;
  obj_value = cplex.getObjValue();
  std::cout << "LR value: " << obj_value << std::endl;

  // Check for integrality
  bool is_integer = true;
  for (int i = 0; i < Xvars.getSize(); ++i)
    if (values[i] < 1 - EPSILON)
      is_integer = false;

  free(mwis_pi);
  COLORstable_freeenv(&mwis_env);
  values.end();
  cplex.end();

  if (is_integer) {
    obj_value = round(obj_value);
    return INTEGER;
  }
  return FRACTIONAL;
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

void LP::branch(std::vector<LP *> &branches, CGVertex cgv) {

  if (N_BRANCHES != 2) {
    std::cout << "N_BRANCHES != 2: Unimplemented" << std::endl;
    abort();
  }

  // Recover hyperedge e and vertex v
  auto [e, v] = cg[cgv];

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

  /* TODO: Remove
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
  */
  return;
}