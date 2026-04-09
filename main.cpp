#include <boost/program_options.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "bp.hpp"
#include "col.hpp"
#include "compact_ilp.hpp"
#include "feas.hpp"
#include "graph.hpp"
#include "heur.hpp"
#include "lp.hpp"
#include "params.hpp"
#include "stats.hpp"

namespace po = boost::program_options;
namespace fs = std::filesystem;

namespace {
struct NullBuffer : std::streambuf {
  int overflow(int c) override { return c; }
};
}  // namespace

// Output directories, initially empty
std::map<std::string, fs::path> outDirs = {
    {"log", fs::path()},    // Directory for log file
    {"debug", fs::path()},  // Directory for debug log file
    {"stat", fs::path()},   // Directory for stats file
    {"sol", fs::path()},    // Directory for solution file
    {"col", fs::path()},    // Directory for column generation log file
    {"iter", fs::path()}    // Directory for iteration log file for semigreedy
};

// Output class to handle std::ostream and file streams
class Output {
 public:
  std::ofstream logFileAux, debugFileAux, colFileAux, iterFileAux;
  std::ostream& logFile;
  std::ostream& debugFile;
  std::ostream& colFile;
  std::ostream& iterFile;

  Output(std::ostream& logStream = std::cout,
         std::ostream& debugStream = std::cout,
         std::ostream& colStream = std::cout,
         std::ostream& iterStream = std::cout)
      : logFileAux(),
        debugFileAux(),
        colFileAux(),
        logFile(logStream),
        debugFile(debugStream),
        colFile(colStream),
        iterFile(iterStream) {}

  Output(std::string logPath, std::string debugPath, std::string colPath,
         std::string iterPath)
      : logFileAux(logPath, std::ofstream::app),
        debugFileAux(debugPath, std::ofstream::app),
        colFileAux(colPath, std::ofstream::app),
        iterFileAux(iterPath, std::ofstream::app),
        logFile(this->logFileAux),
        debugFile(this->debugFileAux),
        colFile(this->colFileAux),
        iterFile(this->iterFileAux) {}
};

int main(int argc, const char** argv) {
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
  desc.add_options()(
      "verbose,v", po::value<int>()->implicit_value(1)->default_value(1),
      "verbosity level (0: quiet, 1: low logging, 2: detailed logging)");
  desc.add_options()("repeat,n", po::value<size_t>()->default_value(1),
                     "number of experiment repetitions");
  desc.add_options()("tree-search", po::value<int>()->default_value(1),
                     "tree search strategy in B&P (1: best-bound, 2: dfs)");
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
  desc.add_options()(
      "heur-2step-variant", po::value<size_t>()->default_value(4),
      "variant of the 2-step heuristic (2: DEG, 3: EDG, 4: AUTO)");
  desc.add_options()("heur-semigreedy-alpha",
                     po::value<double>()->default_value(0.1),
                     "alpha parameter for the semi-greedy heuristic");
  desc.add_options()("heur-semigreedy-iter",
                     po::value<size_t>()->default_value(500),
                     "number of iterations for the semi-greedy heuristic");
  desc.add_options()("feas-root", po::value<int>()->default_value(0),
                     "type of feasibility check for the root node (0: no "
                     "check, 1: enumerative, 2: ILP)");
  desc.add_options()(
      "feas-root-time", po::value<int>()->default_value(300),
      "time limit for feasibility check (in seconds). Only for ILP");
  desc.add_options()("feas-nodes", po::value<int>()->default_value(0),
                     "type of feasibility check for other nodes (0: no "
                     "check, 1: enumerative, 2: ILP)");
  desc.add_options()(
      "feas-nodes-time", po::value<int>()->default_value(60),
      "time limit for feasibility check (in seconds). Only for ILP");
  desc.add_options()(
      "inherit-cols", po::value<int>()->default_value(3),
      "type of column inheritance from parent (0: no inheritance, 1: inherit "
      "all columns, 2: inherit only positive columns, 3: positive columns to "
      "LP, others to pool, 4: all columns to LP)");
  desc.add_options()("dummy-weight", po::value<double>()->default_value(1000.0),
                     "weight of dummy column during initialization");
  desc.add_options()("preproc-off", "do not preprocess the input graph");
  desc.add_options()(
      "pricing-method", po::value<int>()->default_value(6),
      "pricing method (0: ILP, 1: greedy + ILP, 2: greedy + P,Q-MWSSP + ILP, "
      "3: greedy + P-MWSSP + ILP, 4: greedy + P,Q-MWSSP + P-MWSSP + ILP, "
      "5: greedy + P-MWSSP + P,Q-MWSSP + ILP, 6: automatic by density)");
  desc.add_options()("pricing-greedy-max-cols",
                     po::value<size_t>()->default_value(1000),
                     "maximum number of columns to add with greedy pricing");
  desc.add_options()("pricing-max-cols-per-iter",
                     po::value<size_t>()->default_value(10),
                     "maximum number of best columns added per pricing call "
                     "(pool and greedy)");
  desc.add_options()("pricing-greedy-alpha",
                     po::value<double>()->default_value(0.2),
                     "alpha parameter for the greedy pricing heuristic");
  desc.add_options()("pricing-exact-time",
                     po::value<size_t>()->default_value(300),
                     "time limit for exact pricing (in seconds)");
  desc.add_options()("branching-variable", po::value<int>()->default_value(1),
                     "branching variable rule (1: FMS, 2: LNTT)");

  po::positional_options_description pos;
  pos.add("graph", -1);
  po::variables_map vm;

  try {
    po::store(
        po::command_line_parser(argc, argv).options(desc).positional(pos).run(),
        vm);
    po::notify(vm);
  } catch (po::error const& e) {
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
  params.verbose = vm["verbose"].as<int>();
  if (params.verbose < 0 || params.verbose > 2) {
    std::cerr << "verbose must be an integer in [0, 2]" << std::endl;
    return 2;
  }
  params.treeSearch = vm["tree-search"].as<int>();
  if (params.treeSearch < 1 || params.treeSearch > 2) {
    std::cerr << "tree-search must be an integer in [1, 2]" << std::endl;
    return 2;
  }
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
  if (params.inheritColumns < 0 || params.inheritColumns > 4) {
    std::cerr << "inherit-cols must be an integer in [0, 4]" << std::endl;
    return 2;
  }
  params.initializationBigWeight = vm["dummy-weight"].as<double>();
  params.preprocessing = !vm.count("preproc-off");
  params.pricingMethod = vm["pricing-method"].as<int>();
  if (params.pricingMethod < 0 || params.pricingMethod > 6) {
    std::cerr << "pricing-method must be an integer in [0, 6]" << std::endl;
    return 2;
  }
  params.pricingHeur1MaxNCols = vm["pricing-greedy-max-cols"].as<size_t>();
  params.pricingMaxColsPerIter = vm["pricing-max-cols-per-iter"].as<size_t>();
  params.pricingHeur1Alpha = vm["pricing-greedy-alpha"].as<double>();
  params.pricingExactTimeLimit = vm["pricing-exact-time"].as<size_t>();
  params.branchingVariable = vm["branching-variable"].as<int>();
  if (params.branchingVariable < 1 || params.branchingVariable > 2) {
    std::cerr << "branching-variable must be an integer in [1, 2]" << std::endl;
    return 2;
  }

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
  NullBuffer nullBuffer;
  std::ostream nullstream(&nullBuffer);
  for (const std::string& file : inputs) {
    auto path = fs::path(file);
    if (params.is_verbose())
      std::cout << "-> Input file: " << path << std::endl;

    // Open input files
    std::ifstream inGraph(path.string() + ".graph");
    std::ifstream inPartP(path.string() + ".partP");
    std::ifstream inPartQ(path.string() + ".partQ");
    if (!inGraph.is_open() || !inPartP.is_open() || !inPartQ.is_open()) {
      std::cerr << "Error opening files" << std::endl;
      return 2;
    }

    // Read DPCP instance
    Graph graph;
    Partition P, Q;
    read_dpcp_instance(inGraph, inPartP, inPartQ, graph, P, Q);
    size_t nP = P.size();
    size_t nQ = Q.size();

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
      Output out =
          (vm.count("out"))
              ? Output(outDirs["log"].string(), outDirs["debug"].string(),
                       outDirs["col"].string(), outDirs["iter"].string())
              : Output();
      std::ostream& lowLog = params.is_verbose(1) ? out.logFile : nullstream;
      std::ostream& debugLog = params.is_verbose(2) ? out.debugFile : lowLog;

      // Print params
      if (params.is_verbose()) {
        params.print_params(lowLog);
        lowLog << std::endl;
      }

      // Coloring
      Col col;

      // Polymorphic output handling
      auto handle_output = [&](Stats& stats) {
        // Complete stats
        stats.solver = solver;
        stats.instance = path.stem().string();
        stats.run = static_cast<int>(run);
        stats.nvertices = static_cast<int>(num_vertices(graph));
        stats.nedges = static_cast<int>(num_edges(graph));
        stats.nP = static_cast<int>(nP);
        stats.nQ = static_cast<int>(nQ);

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
          col.write_coloring(solFile);
        }
      };

      // Solve the instance
      Stats stats;
      if (solver == "byp") {
        if (params.is_verbose())
          lowLog << "Solving instance " << path << " with B&P" << std::endl;
        DPCPInst dpcp(graph, P, Q);
        BP bp(params, lowLog, debugLog, out.colFile, col,
              vm["ub"].as<double>());
        stats = bp.solve(dpcp);
      } else if (solver == "compact") {
        if (params.is_verbose())
          lowLog << "Solving instance " << path << " with compact ILP"
                 << std::endl;
        DPCPInst dpcp(graph, P, Q);
        stats = solve_ilp(dpcp, params, lowLog, debugLog, col);
      } else if (solver == "heur") {
        if (params.is_verbose())
          lowLog << "Solving instance " << path << " with DPCP heuristic"
                 << std::endl;
        DPCPInst dpcp(graph, P, Q);
        if (params.preprocessing) dpcp.preprocess();
        HeurStats heurStats;
        Col gcol;
        switch (params.heuristicRootNode) {
          case 1:
            heurStats = dpcp_1_step_greedy_heur(dpcp, gcol);
            break;
          case 2:
            heurStats = dpcp_2_step_greedy_heur(dpcp, gcol, params);
            break;
          case 3:
          case 4:
            out.iterFile << path.stem().string() << "," << solver;
            heurStats =
                dpcp_2_step_semigreedy_heur(dpcp, gcol, params, out.iterFile);
            break;
          default:
            std::cerr << "Unknown heuristic root node: "
                      << params.heuristicRootNode << std::endl;
            return 2;
        }
        // Recover coloring for the original graph
        DPCPInst origDpcp(graph, P, Q);
        col = gcol.translate_coloring(dpcp, origDpcp);
        col.color_isolated_vertices(origDpcp, dpcp.get_isolated_vertices());
        assert(col.check_coloring(origDpcp));
        handle_output(heurStats);
        continue;
      } else if (solver == "feas-enum") {
        if (params.is_verbose())
          lowLog << "Deciding feasibility of instance " << path
                 << " with enumerative method" << std::endl;
        DPCPInst dpcp(graph, P, Q);
        dpcp.preprocess();
        stats = dpcp_decide_feasibility_enumerative(dpcp, col, lowLog);
      } else if (solver == "feas-ilp") {
        if (params.is_verbose())
          lowLog << "Deciding feasibility of instance " << path << " with ILP"
                 << std::endl;
        DPCPInst dpcp(graph, P, Q);
        dpcp.preprocess();
        stats =
            dpcp_decide_feasibility_ilp(dpcp, col, params.timeLimit, lowLog);
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