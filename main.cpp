#include "bp.hpp"
#include "col.hpp"
#include "compact_ilp.hpp"
#include "feas.hpp"
#include "graph.hpp"
#include "heur.hpp"
#include "lp.hpp"
#include "params.hpp"
#include "stats.hpp"

#include <boost/program_options.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace po = boost::program_options;
namespace fs = std::filesystem;

// Output directories, initially empty
std::map<std::string, fs::path> outDirs = {
    {"log", fs::path()},  // Directory for log file
    {"stat", fs::path()}, // Directory for stats file
    {"sol", fs::path()},  // Directory for solution file
    {"col", fs::path()}   // Directory for column generation log file
};

// Output class to handle std::ostream and file streams
class Output {
public:
  std::ofstream logFileAux, colFileAux;
  std::ostream &logFile, &colFile;
  Output(std::ostream &logStream = std::cout,
         std::ostream &colStream = std::cout)
      : logFileAux(), colFileAux(), logFile(logStream), colFile(colStream) {}
  Output(std::string logPath, std::string colPath)
      : logFileAux(logPath, std::ofstream::app),
        colFileAux(colPath, std::ofstream::app), logFile(this->logFileAux),
        colFile(this->colFileAux) {}
};

int main(int argc, const char **argv) {

  // Parse arguments
  po::options_description desc(argv[0]);
  desc.add_options()("help,h", "show this help");
  desc.add_options()("solver,s", po::value<std::string>()->required(),
                     "solver, can be any of "
                     "[byp, compact, heur, feas-enum, feas-ilp]");
  desc.add_options()(
      "graph,f",
      po::value<std::vector<std::string>>()->required()->multitoken(),
      "list with the path of each input file (without extension)");
  desc.add_options()(
      "out,o", po::value<std::string>(),
      "output directory (if not given, outputs are printed on stdout)");
  desc.add_options()("time,t", po::value<size_t>()->default_value(900),
                     "time limit in seconds (default: 900)");
  desc.add_options()("repeat,n", po::value<size_t>()->default_value(1),
                     "number of experiment repetitions");
  desc.add_options()("dfs", "use depth-first strategy in the B&P tree");
  desc.add_options()("relax", "solve only the root node");
  desc.add_options()("ub", po::value<double>()->default_value(DBL_MAX),
                     "initial upper bound on the optimal solution value");
  desc.add_options()(
      "heur-root", po::value<int>()->default_value(3),
      "type of heuristic for the root node (0: no heuristic, "
      "1: greedy 1-step, 2: greedy 2-step, 3: semi-greedy 2-step)");
  desc.add_options()("heur-nodes", po::value<int>()->default_value(2),
                     "type of heuristic for other nodes (0: no heuristic, 1: "
                     "greedy 1-step, 2: greedy 2-step, 3: semi-greedy 2-step)");
  desc.add_options()("heur-2step-variant",
                     po::value<size_t>()->default_value(3),
                     "variant of the 2-step heuristic");
  desc.add_options()("heur-semigreedy-alpha",
                     po::value<double>()->default_value(0.1),
                     "alpha parameter for the semi-greedy heuristic");
  desc.add_options()("heur-semigreedy-iter",
                     po::value<size_t>()->default_value(100),
                     "number of iterations for the semi-greedy heuristic");
  desc.add_options()("feas-root", po::value<int>()->default_value(2),
                     "type of feasibility check for the root node (0: no "
                     "check, 1: enumerative, 2: ILP)");
  desc.add_options()(
      "feas-root-time", po::value<int>()->default_value(300),
      "time limit for feasibility check (in seconds). Only for ILP");
  desc.add_options()("feas-nodes", po::value<int>()->default_value(2),
                     "type of feasibility check for other nodes (0: no "
                     "check, 1: enumerative, 2: ILP)");
  desc.add_options()(
      "feas-nodes-time", po::value<int>()->default_value(60),
      "time limit for feasibility check (in seconds). Only for ILP");
  desc.add_options()(
      "inherit-cols", po::value<int>()->default_value(1),
      "type of column inheritance from parent (0: no inheritance, 1: inherit "
      "all columns, 2: inherit only basic columns)");
  desc.add_options()("dummy-weight", po::value<double>()->default_value(1000.0),
                     "weight of dummy column during initialization");
  desc.add_options()("preproc-off", "do not preprocess the input graph");
  desc.add_options()("pool", "use a pool of columns (currently unimplemented)");
  desc.add_options()("pricing-greedy-off",
                     "do not use greedy heuristic for pricing");
  desc.add_options()("pricing-pq-mwsp-off",
                     "do not use P,Q-MWSSP heuristic for pricing");
  desc.add_options()("pricing-p-mwsp-off",
                     "do not use P-MWSSP heuristic for pricing");
  desc.add_options()(
      "pricing-order", po::value<int>()->default_value(1),
      "order of pricing (1: pool -> greedy -> P,Q-MWSSP -> P-MWSSP "
      "II, 2: pool -> greedy -> P-MWSSP II -> P,Q-MWSSP)");
  desc.add_options()("pricing-greedy-max-cols",
                     po::value<size_t>()->default_value(1),
                     "maximum number of columns to add with greedy pricing");
  desc.add_options()("pricing-exact-time",
                     po::value<size_t>()->default_value(300),
                     "time limit for exact pricing (in seconds)");
  desc.add_options()("branching-fms",
                     "use Furini-Malaguti-Santini branching rule");

  po::positional_options_description pos;
  pos.add("graph", -1);
  po::variables_map vm;

  try {
    po::store(
        po::command_line_parser(argc, argv).options(desc).positional(pos).run(),
        vm);
    po::notify(vm);
  } catch (po::error const &e) {
    std::cerr << e.what() << std::endl;
    desc.print(std::cout);
    return 2;
  }

  // Parse parameters
  Params params;
  if (vm.count("help")) {
    desc.print(std::cout);
    return 1;
  }
  std::string solver = vm["solver"].as<std::string>();
  size_t repetitions = vm["repeat"].as<size_t>();
  params.timeLimit = vm["time"].as<size_t>();
  params.dfs = vm.count("dfs");
  params.onlyRelaxation = vm.count("relax");
  params.heuristicRootNode = vm["heur-root"].as<int>();
  params.heuristicOtherNodes = vm["heur-nodes"].as<int>();
  params.heuristic2stepVariant = vm["heur-2step-variant"].as<size_t>();
  params.heuristicSemigreedyAlpha = vm["heur-semigreedy-alpha"].as<double>();
  params.heuristicSemigreedyIter = vm["heur-semigreedy-iter"].as<size_t>();
  params.feasibilityRootNode = vm["feas-root"].as<int>();
  params.feasibilityRootNodeTimeLimit = vm["feas-root-time"].as<int>();
  params.feasibilityOtherNodes = vm["feas-nodes"].as<int>();
  params.feasibilityOtherNodesTimeLimit = vm["feas-nodes-time"].as<int>();
  params.inheritColumns = vm["inherit-cols"].as<int>();
  params.initializationBigWeight = vm["dummy-weight"].as<double>();
  params.preprocStep1 = !vm.count("preproc-off");
  params.preprocStep2 = !vm.count("preproc-off");
  params.preprocStep3 = !vm.count("preproc-off");
  params.preprocStep4 = !vm.count("preproc-off");
  params.usePool = vm.count("pool");
  params.pricingHeur1 = !vm.count("pricing-greedy-off");
  params.pricingHeur2 = !vm.count("pricing-pq-mwsp-off");
  params.pricingHeur3 = !vm.count("pricing-p-mwsp-off");
  params.pricingOrder = vm["pricing-order"].as<int>();
  params.pricingHeur1MaxNCols = vm["pricing-greedy-max-cols"].as<size_t>();
  params.pricingExactTimeLimit = vm["pricing-exact-time"].as<size_t>();
  params.branchingFMS = vm.count("branching-fms");

  // Parse input files
  std::vector<std::string> inputs;
  if (vm.count("graph")) {
    inputs = vm["graph"].as<std::vector<std::string>>();
  } else {
    std::cerr << "no input\n" << std::endl;
    desc.print(std::cout);
    return 2;
  }

  // Parse output directory
  fs::path outDir;
  if (vm.count("out")) {
    outDir = fs::path(vm["out"].as<std::string>());
    // Create output directories
    if (!fs::is_directory(outDir)) {
      fs::create_directory(outDir);
    }
    for (auto [key, _] : outDirs) {
      outDirs[key] = outDir / key / "";
      if (!fs::is_directory(outDirs[key])) {
        fs::create_directories(outDirs[key]);
      }
    }
  }

  // Read inputs
  for (const std::string &file : inputs) {

    auto path = fs::path(file);
    std::cout << "-> Input file: " << path << std::endl;

    // Open input files
    std::ifstream inGraph(path.string() + ".graph");
    std::ifstream inPartA(path.string() + ".partA");
    std::ifstream inPartB(path.string() + ".partB");
    if (!inGraph.is_open() || !inPartA.is_open() || !inPartB.is_open()) {
      std::cerr << "Error opening files" << std::endl;
      return 2;
    }

    // Read DPCP instance
    Graph graph;
    size_t nA, nB;
    std::tie(graph, nA, nB) = read_dpcp_instance(inGraph, inPartA, inPartB);

    for (size_t run = 0; run < repetitions; ++run) {

      // Get name for the current run
      fs::path currentName(path.stem().string() + "-" + solver + "-" +
                           std::to_string(run) + ".txt");
      if (vm.count("out")) {

        // Set output file names
        for (auto [key, _] : outDirs) {
          outDirs[key].replace_filename(currentName);
          outDirs[key].replace_extension(key);
        }

        // Resume check
        if (fs::is_regular_file(outDirs["stat"])) {
          std::ifstream stats;
          stats.open(outDirs["stat"]);
          std::string line;
          while (getline(stats, line)) {
            std::cout << line << std::endl;
          }
          continue;
        }
      }

      // Prepare log and col output files
      Output out = (vm.count("out")) ? Output(outDirs["log"].string(),
                                              outDirs["col"].string())
                                     : Output();

      // Print params
      params.print_params(out.logFile);
      out.logFile << std::endl;

      // Coloring
      Col col;

      // Polymorphic output handling
      auto handle_output = [&](Stats &stats) {
        // Complete stats
        stats.solver = solver;
        stats.instance = path.stem().string();
        stats.run = static_cast<int>(run);
        stats.nvertices = static_cast<int>(num_vertices(graph));
        stats.nedges = static_cast<int>(num_edges(graph));
        stats.nA = static_cast<int>(nA);
        stats.nB = static_cast<int>(nB);

        // Write stats
        if (vm.count("out")) {
          std::ofstream statFile(outDirs["stat"].string(), std::ofstream::app);
          stats.write_stats(statFile);
          stats.print_stats(out.logFile);
        } else
          stats.print_stats(std::cout);
        stats.write_stats(std::cout);
        out.logFile << std::endl;

        // Write coloring
        if (stats.state == OPTIMAL || stats.state == FEASIBLE) {
          std::ofstream solFile(outDirs["sol"].string(), std::ofstream::app);
          col.write_coloring(graph, solFile);
        }
      };

      // Solve the instance
      Stats stats;
      if (solver == "byp") {
        out.logFile << "Solving instance " << path << " with B&P" << std::endl;
        Graph *gcopy = new Graph;
        graph_copy(graph, *gcopy);
        Pool pool;
        LP *lp = new LP(gcopy, params, pool, graph, out.colFile, true);
        Node *root = new Node(lp);
        BP<Col> bp(params, out.logFile, col, vm["ub"].as<double>());
        stats = bp.solve(root);
      } else if (solver == "compact") {
        out.logFile << "Solving instance " << path << " with compact ILP"
                    << std::endl;
        stats = solve_ilp(graph, params, out.logFile, col);
      } else if (solver == "heur") {
        out.logFile << "Solving instance " << path << " with DPCP heuristic"
                    << std::endl;
        GraphEnv genv(&graph, params.preprocStep1, params.preprocStep2,
                      params.preprocStep3, params.preprocStep4, false);
        HeurStats heurStats;
        switch (params.heuristicRootNode) {
        case 1:
          heurStats = dpcp_1_step_greedy_heur(genv, col);
          break;
        case 2:
          heurStats = dpcp_2_step_greedy_heur(genv, col, params);
          break;
        case 3:
          heurStats = dpcp_2_step_semigreedy_heur(genv, col, params);
          break;
        default:
          std::cerr << "Unknown heuristic root node: "
                    << params.heuristicRootNode << std::endl;
          return 2;
        }
        handle_output(heurStats);
        return 0;
      } else if (solver == "feas-enum") {
        out.logFile << "Deciding feasibility of instance " << path
                    << " with enumerative method" << std::endl;
        stats = dpcp_decide_feasibility_enumerative(
            GraphEnv(&graph, false, false, false, false, true), col,
            out.logFile);
      } else if (solver == "feas-ilp") {
        out.logFile << "Deciding feasibility of instance " << path
                    << " with ILP" << std::endl;
        stats = dpcp_decide_feasibility_ilp(
            GraphEnv(&graph, false, false, false, false, true), col,
            params.timeLimit, out.logFile);
      } else {
        std::cerr << "Unknown solver: " << solver << std::endl;
        return 2;
      }

      // Handle output
      handle_output(stats);
    }
  }
  return 0;
}